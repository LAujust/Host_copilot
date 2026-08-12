import math

import pytest

from host_copilot.models import (
    CircleLocalization,
    EllipseLocalization,
    SearchConfig,
    TransientContext,
)


def test_circle_uses_90_percent_containment_by_default():
    localization = CircleLocalization(10.0, 20.0, 10.0)
    assert localization.confidence == 0.9
    assert localization.sigma_arcsec == pytest.approx(
        10.0 / math.sqrt(-2.0 * math.log(0.1))
    )
    assert localization.contains(10.0, 20.0)
    assert localization.relative_likelihood(10.0, 20.0) == pytest.approx(1.0)


def test_circle_handles_ra_wrap():
    localization = CircleLocalization(359.999, 0.0, 10.0)
    assert localization.contains(0.001, 0.0)


def test_ellipse_contains_center_and_rejects_outside_major_axis():
    localization = EllipseLocalization(120.0, -30.0, 10.0, 5.0, 0.0)
    assert localization.contains(120.0, -30.0)
    # Approximately 12 arcsec north, along PA=0 major axis.
    assert not localization.contains(120.0, -30.0 + 12.0 / 3600.0)


def test_transient_requires_complete_optical_position():
    with pytest.raises(ValueError, match="supplied together"):
        TransientContext(CircleLocalization(10.0, 20.0, 10.0), optical_ra_deg=10.0)


def test_mode_defaults_are_distinct():
    quick = SearchConfig(mode="quick")
    full = SearchConfig(mode="full")
    assert quick.z_max == 0.1
    assert quick.deadline_seconds == 30.0
    assert quick.providers == ("regalade",)
    assert full.z_max == 0.5
    assert full.deadline_seconds == 180.0
    assert "legacy" in full.providers
    assert full.image_recovery is True


def test_explicit_query_geometry_overrides_localization_margin():
    transient = TransientContext(CircleLocalization(10.0, 20.0, 10.0))
    config = SearchConfig(
        mode="full",
        search_ra_deg=10.1,
        search_dec_deg=20.1,
        search_radius_arcsec=20.0,
    )
    assert config.query_geometry(transient) == pytest.approx((10.1, 20.1, 20.0))


def test_query_geometry_override_must_be_complete():
    with pytest.raises(ValueError, match="must be supplied together"):
        SearchConfig(mode="full", search_ra_deg=10.0)
