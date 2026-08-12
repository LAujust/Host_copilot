"""Optional image-based recovery of catalog-missed host candidates."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any
import warnings

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import SigmaClip
from astropy.wcs import FITSFixedWarning, WCS
from astropy.wcs.utils import proj_plane_pixel_scales
import numpy as np
import pandas as pd
import requests

from .association import stable_candidate_id
from .image import Imager
from .models import HostCandidate, SearchConfig, TransientContext


class ImageRecovery:
    """Detect extended sources in a survey cutout with Photutils.

    This is a recovery path, not a replacement for survey model catalogs.  Its
    candidates deliberately carry no redshift and are ranked with the
    full-mode unknown-redshift prior.
    """

    def __init__(self, save_path: str | Path):
        self.save_path = Path(save_path)
        self.save_path.mkdir(parents=True, exist_ok=True)

    def _cutout_path(
        self, transient: TransientContext, config: SearchConfig
    ) -> Path | None:
        ra, dec, query_radius = config.query_geometry(transient)
        for filename in ("ps1_r_ref.fits", "ls_r.fits", "ls_r_science.fits"):
            candidate = self.save_path / filename
            if candidate.exists() and self._covers_search(
                candidate, ra, dec, query_radius
            ):
                return candidate
        loc = transient.localization
        if config.search_radius_arcsec is not None:
            ra, dec, radius = config.query_geometry(transient)
        else:
            ra, dec = loc.ra_deg, loc.dec_deg
            radius = loc.enclosing_radius_arcsec + min(
                config.association_margin_arcsec, 60.0
            )
        size_pixels = max(64, int(math.ceil(2.0 * radius / 0.262)))

        def download_legacy(product: str, layer: str) -> Path:
            response = requests.get(
                "https://www.legacysurvey.org/viewer/fits-cutout",
                params={
                    "ra": ra,
                    "dec": dec,
                    "layer": layer,
                    "pixscale": 0.262,
                    "bands": "r",
                    "size": size_pixels,
                },
                timeout=float(config.provider_timeout_seconds),
            )
            response.raise_for_status()
            product_path = self.save_path / f"ls_r_{product}.fits"
            product_path.write_bytes(response.content)
            self._read_image(product_path)
            return product_path

        try:
            science_path = download_legacy("science", "ls-dr10")
            for product, layer in (
                ("model", "ls-dr10-model"),
                ("residual", "ls-dr10-resid"),
            ):
                try:
                    download_legacy(product, layer)
                except Exception:
                    continue
            return science_path
        except Exception:
            pass

        # Pan-STARRS is the northern-sky fallback.  Call its cutout method
        # directly so a second Legacy Survey attempt is not made.
        try:
            path = Imager(
                ra,
                dec,
                radius,
                band="r",
                save_path=str(self.save_path),
            ).PS_cutout()
            return None if path is None else Path(path)
        except Exception:
            return None

    def get_reference_cutout(
        self, transient: TransientContext, config: SearchConfig
    ) -> Path | None:
        """Return a cached or newly downloaded science reference FITS cutout."""

        return self._cutout_path(transient, config)

    @classmethod
    def _covers_search(
        cls, path: Path, ra_deg: float, dec_deg: float, radius_arcsec: float
    ) -> bool:
        """Check that a cached cutout contains the requested cone."""

        try:
            data, header = cls._read_image(path)
            wcs = cls._celestial_wcs(header)
            x, y = wcs.world_to_pixel(SkyCoord(ra_deg * u.deg, dec_deg * u.deg))
            scales = proj_plane_pixel_scales(wcs) * u.deg
            radius_pixels = radius_arcsec / float(np.min(scales.to_value(u.arcsec)))
            height, width = data.shape
            return (
                radius_pixels <= x + 0.5 < width - radius_pixels + 0.5
                and radius_pixels <= y + 0.5 < height - radius_pixels + 0.5
            )
        except Exception:
            return False

    @staticmethod
    def _read_image(path: Path) -> tuple[np.ndarray, fits.Header]:
        with fits.open(path, memmap=False) as hdul:
            for hdu in hdul:
                if hdu.data is None:
                    continue
                data = np.asarray(hdu.data, dtype=float)
                while data.ndim > 2:
                    data = data[0]
                if data.ndim == 2:
                    return data, hdu.header.copy()
        raise ValueError(f"no two-dimensional image found in {path}")

    @staticmethod
    def _celestial_wcs(header: fits.Header) -> WCS:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            return WCS(header).celestial

    def run(
        self, transient: TransientContext, config: SearchConfig
    ) -> tuple[list[HostCandidate], dict[str, Any]]:
        try:
            from photutils.background import Background2D, MedianBackground
            from photutils.segmentation import (
                SourceCatalog,
                deblend_sources,
                detect_sources,
                detect_threshold,
            )
        except ImportError as exc:
            return [], {"state": "unavailable", "error": str(exc)}

        try:
            path = self._cutout_path(transient, config)
            if path is None:
                return [], {"state": "no_cutout"}
            data, header = self._read_image(path)
            wcs = self._celestial_wcs(header)
            mask = ~np.isfinite(data)
            finite_pixels = data[~mask]
            if finite_pixels.size < 100:
                raise ValueError("cutout has too few finite pixels")

            minimum_dimension = min(data.shape)
            box = max(16, min(64, minimum_dimension // 4))
            background = Background2D(
                data,
                (box, box),
                filter_size=(3, 3),
                sigma_clip=SigmaClip(sigma=3.0),
                bkg_estimator=MedianBackground(),
                mask=mask,
            )
            data_sub = data - background.background
            threshold = detect_threshold(
                data_sub,
                nsigma=2.5,
                background=0.0,
                error=background.background_rms,
                mask=mask,
            )
            segmentation = detect_sources(data_sub, threshold, npixels=8, mask=mask)
            if segmentation is None:
                return [], {"state": "empty", "cutout": str(path)}
            try:
                segmentation = deblend_sources(
                    data_sub,
                    segmentation,
                    npixels=8,
                    nlevels=32,
                    contrast=0.001,
                    progress_bar=False,
                )
            except Exception:
                # A valid non-deblended segmentation is preferable to losing
                # the recovery result entirely.
                pass

            catalog = SourceCatalog(
                data_sub,
                segmentation,
                background=background.background,
                error=background.background_rms,
                mask=mask,
                wcs=wcs,
            )
            table = catalog.to_table(
                columns=(
                    "label",
                    "xcentroid",
                    "ycentroid",
                    "sky_centroid",
                    "area",
                    "semimajor_sigma",
                    "semiminor_sigma",
                    "orientation",
                    "segment_flux",
                )
            )
            scales = proj_plane_pixel_scales(wcs) * u.deg
            pixel_scale_arcsec = float(np.mean(scales.to_value(u.arcsec)))
            candidates: list[HostCandidate] = []
            records: list[dict[str, Any]] = []
            for row in table:
                sky = row["sky_centroid"]
                if sky is None:
                    continue
                area_value = float(getattr(row["area"], "value", row["area"]))
                if area_value < 8:
                    continue
                major_pix = float(
                    getattr(row["semimajor_sigma"], "value", row["semimajor_sigma"])
                )
                minor_pix = float(
                    getattr(row["semiminor_sigma"], "value", row["semiminor_sigma"])
                )
                # Three Gaussian sigma encloses most of a smooth detection.
                major = max(pixel_scale_arcsec, 3.0 * major_pix * pixel_scale_arcsec)
                minor = max(pixel_scale_arcsec, 3.0 * minor_pix * pixel_scale_arcsec)
                orientation = float(
                    getattr(row["orientation"], "value", row["orientation"])
                )
                label = int(row["label"])
                source_id = f"{path.name}:{label}"
                candidate = HostCandidate(
                    candidate_id=stable_candidate_id(
                        "image", source_id, float(sky.ra.deg), float(sky.dec.deg)
                    ),
                    ra_deg=float(sky.ra.deg),
                    dec_deg=float(sky.dec.deg),
                    name=f"image source {label}",
                    catalogs=["image"],
                    catalog_ids={"image": source_id},
                    semimajor_arcsec=major,
                    semiminor_arcsec=minor,
                    position_angle_deg=orientation,
                    morphology="image_only",
                    quality_flags=["image_only", "redshift_unknown"],
                    extra={
                        "segment_area_pix2": area_value,
                        "segment_flux": float(row["segment_flux"]),
                        "cutout": str(path),
                    },
                )
                candidates.append(candidate)
                records.append(candidate.to_record())

            segmentation_path = self.save_path / "segmentation.fits"
            fits.PrimaryHDU(
                np.asarray(segmentation.data, dtype=np.int32), header=header
            ).writeto(segmentation_path, overwrite=True)
            pd.DataFrame(records).to_csv(
                self.save_path / "image_candidates.csv", index=False
            )
            overlay_path: Path | None = self.save_path / "image_detection_overlay.png"
            try:
                import matplotlib.pyplot as plt

                figure, axis = plt.subplots(figsize=(8, 8))
                low, high = np.nanpercentile(data, (1, 99))
                axis.imshow(data, origin="lower", cmap="gray", vmin=low, vmax=high)
                for row in table:
                    axis.add_patch(
                        plt.Circle(
                            (float(row["xcentroid"]), float(row["ycentroid"])),
                            radius=6,
                            fill=False,
                            edgecolor="cyan",
                            linewidth=0.8,
                        )
                    )
                axis.set_title("Photutils source-recovery detections")
                figure.tight_layout()
                figure.savefig(overlay_path, dpi=150)
                plt.close(figure)
            except Exception:
                overlay_path = None
            return candidates, {
                "state": "success" if candidates else "empty",
                "cutout": str(path),
                "segmentation": str(segmentation_path),
                "overlay": str(overlay_path) if overlay_path is not None else None,
                "row_count": len(candidates),
            }
        except Exception as exc:
            return [], {"state": "error", "error": f"{type(exc).__name__}: {exc}"}
