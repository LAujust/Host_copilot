"""High-level quick and full host-galaxy search pipeline."""

from __future__ import annotations

from concurrent.futures import Future, ThreadPoolExecutor, as_completed
import math
from pathlib import Path
from threading import BoundedSemaphore
import time
from typing import Any

import astropy.units as u
from astropy.coordinates import Angle
from astropy.table import Table

from .association import (
    apply_gaia_star_veto,
    deduplicate_candidates,
    filter_by_redshift,
    rank_candidates,
)
from .cache import QueryCache
from .catalog import GalaxyFinder
from .extraction import ImageRecovery
from .image import Imager
from .models import (
    CircleLocalization,
    EllipseLocalization,
    HostCandidate,
    HostSearchResult,
    ProviderStatus,
    SearchConfig,
    TransientContext,
)
from .providers import CatalogProvider, ProviderResult, build_providers


class HostPipeline:
    """Search for and rank host-galaxy candidates.

    The positional constructor is retained for compatibility.  New code should
    prefer ``localization=...`` and :meth:`search`, which returns a structured
    :class:`HostSearchResult`.
    """

    def __init__(
        self,
        ra: float | None = None,
        dec: float | None = None,
        r_arcsec: float | None = None,
        zcutout: float | None = None,
        quick: bool = True,
        save_path: str = "./",
        *,
        localization: CircleLocalization | EllipseLocalization | None = None,
        mode: str | None = None,
        optical_ra: float | None = None,
        optical_dec: float | None = None,
        transient_redshift: float | None = None,
        transient_redshift_error: float | None = None,
        transient_name: str = "",
        classification: str | None = None,
        cache_dir: str | Path | None = None,
    ):
        if localization is None:
            if ra is None or dec is None or r_arcsec is None:
                raise ValueError("provide localization or all of ra, dec, and r_arcsec")
            localization = CircleLocalization(float(ra), float(dec), float(r_arcsec))
        self.localization = localization
        self.ra = localization.ra_deg
        self.dec = localization.dec_deg
        self.r_arcsec = localization.enclosing_radius_arcsec
        self.quick = quick if mode is None else mode == "quick"
        self.mode = mode or ("quick" if quick else "full")
        if self.mode not in {"quick", "full"}:
            raise ValueError("mode must be 'quick' or 'full'")
        self.zcutout = zcutout if zcutout is not None else (0.1 if self.quick else 0.5)
        self.save_path = Path(save_path)
        self.save_path.mkdir(parents=True, exist_ok=True)
        self.cache_dir = (
            Path(cache_dir) if cache_dir else self.save_path / ".host_copilot_cache"
        )
        self.context = TransientContext(
            localization=localization,
            name=transient_name,
            optical_ra_deg=optical_ra,
            optical_dec_deg=optical_dec,
            redshift=transient_redshift,
            redshift_error=transient_redshift_error,
            classification=classification,
        )
        self.pos = localization.center

        # Compatibility attributes used by existing notebooks and batch code.
        self.galaxy_finder = GalaxyFinder(
            self.ra, self.dec, self.r_arcsec, save_path=str(self.save_path)
        )
        self.imager = Imager(
            self.ra, self.dec, self.r_arcsec, save_path=str(self.save_path)
        )
        self.last_result: HostSearchResult | None = None
        self._datalab_slots = BoundedSemaphore(2)

    def _query_provider(
        self,
        provider: CatalogProvider,
        transient: TransientContext,
        config: SearchConfig,
        cache: QueryCache,
    ) -> ProviderResult:
        if getattr(provider, "tap_url", "").startswith("https://datalab.noirlab.edu"):
            with self._datalab_slots:
                return provider.query(transient, config, cache)
        return provider.query(transient, config, cache)

    def _run_providers(
        self,
        transient: TransientContext,
        config: SearchConfig,
        cache: QueryCache,
    ) -> dict[str, ProviderResult]:
        providers = build_providers(config.providers or ())
        if config.mode == "quick" or len(providers) <= 1:
            return {
                provider.name: self._query_provider(provider, transient, config, cache)
                for provider in providers
            }

        results: dict[str, ProviderResult] = {}
        executor = ThreadPoolExecutor(max_workers=min(6, len(providers)))
        futures: dict[Future[ProviderResult], CatalogProvider] = {
            executor.submit(
                self._query_provider, provider, transient, config, cache
            ): provider
            for provider in providers
        }
        try:
            for future in as_completed(futures, timeout=config.deadline_seconds):
                provider = futures[future]
                try:
                    results[provider.name] = future.result()
                except (
                    Exception
                ) as exc:  # Defensive; providers normally contain errors.
                    results[provider.name] = ProviderResult(
                        [],
                        ProviderStatus(
                            provider=provider.name,
                            state="service_error",
                            error=f"{type(exc).__name__}: {exc}",
                            catalog_version=provider.version,
                        ),
                    )
        except TimeoutError:
            pass
        finally:
            for future, provider in futures.items():
                if provider.name in results:
                    continue
                future.cancel()
                results[provider.name] = ProviderResult(
                    [],
                    ProviderStatus(
                        provider=provider.name,
                        state="timeout",
                        elapsed_seconds=float(config.deadline_seconds),
                        catalog_version=provider.version,
                        error="global full-mode deadline exceeded",
                    ),
                )
            executor.shutdown(wait=False, cancel_futures=True)
        return results

    def search(
        self,
        transient: TransientContext | None = None,
        config: SearchConfig | None = None,
    ) -> HostSearchResult:
        """Run a structured quick or full host search."""

        started = time.monotonic()
        transient = transient or self.context
        config = config or SearchConfig(
            mode=self.mode,
            z_max=self.zcutout,
            cache_dir=self.cache_dir,
        )
        if config.cache_dir is None:
            config.cache_dir = self.cache_dir
        cache = QueryCache(config.cache_dir)
        provider_results = self._run_providers(transient, config, cache)

        statuses: dict[str, ProviderStatus] = {}
        host_candidates: list[HostCandidate] = []
        gaia_sources: list[HostCandidate] = []
        warnings: list[str] = []
        for provider_name, provider_result in provider_results.items():
            if provider_result.status is not None:
                statuses[provider_name] = provider_result.status
                if provider_result.status.warning:
                    warnings.append(
                        f"{provider_name}: {provider_result.status.warning}"
                    )
                if provider_result.status.state not in {
                    "success",
                    "empty",
                    "cached",
                }:
                    warnings.append(f"{provider_name}: {provider_result.status.state}")
            if provider_name == "gaia":
                gaia_sources.extend(provider_result.candidates)
            else:
                host_candidates.extend(provider_result.candidates)

        if config.image_recovery:
            image_started = time.monotonic()
            image_candidates, image_info = ImageRecovery(self.save_path).run(
                transient, config
            )
            host_candidates.extend(image_candidates)
            image_state = str(image_info.get("state", "error"))
            status_state = (
                "success"
                if image_state == "success"
                else "empty"
                if image_state in {"empty", "no_cutout"}
                else "service_error"
            )
            statuses["image"] = ProviderStatus(
                provider="image",
                state=status_state,
                elapsed_seconds=time.monotonic() - image_started,
                row_count=len(image_candidates),
                warning=(
                    None
                    if status_state in {"success", "empty"}
                    else str(image_info.get("error", image_state))
                ),
            )
            if status_state == "service_error":
                warnings.append(f"image: {image_info.get('error', image_state)}")
        else:
            image_info = {"state": "disabled"}

        host_candidates = deduplicate_candidates(host_candidates)
        apply_gaia_star_veto(host_candidates, gaia_sources)
        host_candidates = filter_by_redshift(host_candidates, config)
        host_candidates = rank_candidates(host_candidates, transient, config)

        result = HostSearchResult(
            transient=transient,
            config=config,
            candidates=host_candidates,
            provider_status=statuses,
            warnings=warnings,
            metadata={
                "elapsed_seconds": time.monotonic() - started,
                "ranking_calibrated": False,
                "relative_probability_note": (
                    "relative_probability is an uncalibrated normalized score; "
                    "posterior_probability is intentionally unset"
                ),
                "image_recovery": image_info,
                "catalog_candidate_count_before_deduplication": sum(
                    len(item.candidates)
                    for name, item in provider_results.items()
                    if name != "gaia"
                ),
                "gaia_reference_count": len(gaia_sources),
            },
        )
        if not result.complete:
            result.warnings.append(
                "Search is partial or truncated; do not claim catalog completeness."
            )
        self.last_result = result
        self._update_compatibility_catalog(result)
        return result

    def _update_compatibility_catalog(self, result: HostSearchResult) -> None:
        records: list[dict[str, Any]] = []
        for candidate in result.candidates:
            if "regalade" not in candidate.catalogs:
                continue
            records.append(
                {
                    "Name": candidate.name,
                    "RAJ2000": candidate.ra_deg,
                    "DEJ2000": candidate.dec_deg,
                    "z": candidate.redshift,
                    "R1": candidate.semimajor_arcsec,
                    "R2": candidate.semiminor_arcsec,
                    "PA": candidate.position_angle_deg,
                    "sep": candidate.separation_arcsec,
                }
            )
        import pandas as pd

        self.galaxy_finder.reglade_df = pd.DataFrame(records)

    def _candidate_visualization_table(self, result: HostSearchResult) -> Table:
        rows = []
        for candidate in result.candidates:
            rows.append(
                (
                    candidate.name,
                    candidate.ra_deg,
                    candidate.dec_deg,
                    candidate.redshift if candidate.redshift is not None else math.nan,
                    candidate.semimajor_arcsec or 1.0,
                    candidate.semiminor_arcsec or 1.0,
                    candidate.position_angle_deg or 0.0,
                    candidate.separation_arcsec or 0.0,
                    candidate.rank or 0,
                    candidate.relative_probability,
                    candidate.assessment,
                )
            )
        table = Table(
            rows=rows,
            names=(
                "Name",
                "RAJ2000",
                "DEJ2000",
                "z",
                "R1",
                "R2",
                "PA",
                "sep",
                "rank",
                "relative_probability",
                "assessment",
            ),
        )
        table["R1"].unit = u.arcsec
        table["R2"].unit = u.arcsec
        table["PA"].unit = u.deg
        table["sep"].unit = u.arcsec
        return table

    def build_aladin(self, result: HostSearchResult | None = None) -> Any:
        """Build an ipyaladin visualization for a structured result."""

        from ipyaladin import Aladin, EllipseError
        from regions import CircleSkyRegion, EllipseSkyRegion

        result = result or self.last_result
        if result is None:
            raise RuntimeError("search must run before visualization")
        loc = result.transient.localization
        fov = max(0.02, 3.0 * loc.enclosing_radius_arcsec / 3600.0)
        aladin = Aladin(
            fov=fov,
            target=loc.center,
            survey="CDS/P/PanSTARRS/DR1/color-z-zg-g",
        )
        table = self._candidate_visualization_table(result)
        if len(table):
            aladin.add_table(
                table,
                shape=EllipseError(
                    maj_axis="R1",
                    min_axis="R2",
                    angle="PA",
                    default_shape="cross",
                ),
                color="cyan",
            )
        if isinstance(loc, CircleLocalization):
            region = CircleSkyRegion(
                center=loc.center,
                radius=Angle(loc.radius_arcsec, "arcsec"),
                visual={"edgecolor": "yellow", "linestyle": "dashed"},
            )
        else:
            region = EllipseSkyRegion(
                center=loc.center,
                width=Angle(2.0 * loc.semimajor_arcsec, "arcsec"),
                height=Angle(2.0 * loc.semiminor_arcsec, "arcsec"),
                angle=Angle(loc.position_angle_deg, "deg"),
                visual={"edgecolor": "yellow", "linestyle": "dashed"},
            )
        aladin.add_graphic_overlay_from_region([region])
        result.aladin = aladin
        return aladin

    def filter_and_visualize(self) -> tuple[Any, Table | None]:
        """Compatibility method returning the historical tuple."""

        result = self.last_result or self.search()
        aladin = self.build_aladin(result)
        table = self._candidate_visualization_table(result)
        return aladin, table if len(table) else None

    def run(self) -> tuple[Any, Table | None]:
        """Run the selected mode and return the historical ``(Aladin, Table)``."""

        print("=" * 50)
        print(f"[{self.mode.upper()} MODE]")
        print("Searching for galaxies...")
        result = self.search()
        for candidate in result.candidates:
            z_text = (
                "unknown" if candidate.redshift is None else f"{candidate.redshift:.4f}"
            )
            sep_text = (
                "unknown"
                if candidate.separation_arcsec is None
                else f'{candidate.separation_arcsec:.2f}"'
            )
            print(
                f"#{candidate.rank} {candidate.name or candidate.candidate_id}: "
                f"z={z_text}, sep={sep_text}, {candidate.assessment}"
            )
        if result.warnings:
            print("Warnings:")
            for warning in result.warnings:
                print(f"- {warning}")
        aladin = self.build_aladin(result)
        table = self._candidate_visualization_table(result)
        print("HostPipeline run completed.")
        print("=" * 50)
        return aladin, table if len(table) else None
