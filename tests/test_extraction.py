import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from host_copilot.extraction import ImageRecovery
from host_copilot.models import CircleLocalization, SearchConfig, TransientContext


def test_image_recovery_creates_unknown_redshift_candidates(tmp_path):
    rng = np.random.default_rng(1234)
    yy, xx = np.mgrid[:128, :128]
    image = rng.normal(100.0, 1.0, size=(128, 128))
    image += 80.0 * np.exp(-0.5 * (((xx - 64.0) / 5.0) ** 2 + ((yy - 64.0) / 3.0) ** 2))

    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [64.0, 64.0]
    wcs.wcs.cdelt = np.array([-0.262 / 3600.0, 0.262 / 3600.0])
    wcs.wcs.crval = [10.0, 20.0]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    fits.PrimaryHDU(image, header=wcs.to_header()).writeto(tmp_path / "ps1_r_ref.fits")

    candidates, status = ImageRecovery(tmp_path).run(
        TransientContext(CircleLocalization(10.0, 20.0, 10.0)),
        SearchConfig(
            mode="full",
            image_recovery=True,
            search_ra_deg=10.0,
            search_dec_deg=20.0,
            search_radius_arcsec=10.0,
        ),
    )
    assert status["state"] == "success"
    assert candidates
    assert all(candidate.redshift is None for candidate in candidates)
    assert any("image_only" in candidate.quality_flags for candidate in candidates)
    assert (tmp_path / "segmentation.fits").exists()
    assert (tmp_path / "image_candidates.csv").exists()


def test_cached_cutout_must_cover_new_search_center(tmp_path):
    image = np.ones((100, 100))
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [50.0, 50.0]
    wcs.wcs.cdelt = np.array([-1.0 / 3600.0, 1.0 / 3600.0])
    wcs.wcs.crval = [10.0, 20.0]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    path = tmp_path / "cutout.fits"
    fits.PrimaryHDU(image, header=wcs.to_header()).writeto(path)

    assert ImageRecovery._covers_search(path, 10.0, 20.0, 20.0)
    assert not ImageRecovery._covers_search(path, 10.02, 20.0, 20.0)
