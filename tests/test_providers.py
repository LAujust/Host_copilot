from types import SimpleNamespace

from host_copilot.cache import QueryCache
from host_copilot.models import CircleLocalization, SearchConfig, TransientContext
from host_copilot.providers import (
    DesiProvider,
    GaiaProvider,
    LegacySurveyProvider,
    NedProvider,
    PanStarrsProvider,
    RegaladeProvider,
)


def test_regalade_normalization_preserves_shape_and_redshift_flags():
    rows = [
        {
            "Name": "Galaxy A",
            "RAJ2000": 10.0,
            "DEJ2000": 20.0,
            "z": 0.05,
            "z_spec": 1,
            "R1": 8.0,
            "R2": 4.0,
            "PA": 45.0,
            "fRel": 1,
        }
    ]
    result = RegaladeProvider().normalize(rows)
    assert len(result) == 1
    assert result[0].redshift_kind == "spec"
    assert result[0].semimajor_arcsec == 8.0
    assert "low_reliability" in result[0].quality_flags


def test_regalade_comment_only_response_is_empty_not_error(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "host_copilot.providers._bounded_request",
        lambda *args, **kwargs: SimpleNamespace(text="# VizieR\n# no rows\n"),
    )
    context = TransientContext(CircleLocalization(10.0, 20.0, 10.0))
    config = SearchConfig(mode="quick", image_recovery=False, cache_dir=tmp_path)
    result = RegaladeProvider().query(context, config, QueryCache(tmp_path))
    assert result.candidates == []
    assert result.status is not None
    assert result.status.state == "empty"
    assert result.status.error is None


def test_panstarrs_blank_response_is_empty(monkeypatch):
    monkeypatch.setattr(
        "host_copilot.providers._bounded_request",
        lambda *args, **kwargs: SimpleNamespace(text="\n"),
    )
    context = TransientContext(CircleLocalization(10.0, 20.0, 10.0))
    rows = PanStarrsProvider().fetch(
        context, SearchConfig(mode="full", image_recovery=False)
    )
    assert rows == []


def test_tap_header_only_response_is_empty(monkeypatch):
    monkeypatch.setattr(
        "host_copilot.providers._bounded_request",
        lambda *args, **kwargs: SimpleNamespace(text="ls_id,ra,dec\n"),
    )
    context = TransientContext(CircleLocalization(10.0, 20.0, 10.0))
    rows = LegacySurveyProvider().fetch(
        context, SearchConfig(mode="full", image_recovery=False)
    )
    assert rows == []


def test_legacy_photo_z_and_shape_normalization():
    rows = [
        {
            "ls_id": 123,
            "ra": 10.0,
            "dec": 20.0,
            "type": "EXP",
            "shape_r": 2.0,
            "shape_e1": 0.2,
            "shape_e2": 0.0,
            "z_spec": None,
            "z_phot_median": 0.3,
            "z_phot_std": 0.05,
            "mag_r": 21.0,
        }
    ]
    result = LegacySurveyProvider().normalize(rows)
    assert result[0].redshift_kind == "photo"
    assert result[0].redshift_error == 0.05
    assert result[0].semiminor_arcsec < result[0].semimajor_arcsec


def test_gaia_requires_significant_astrometry_or_high_star_probability():
    rows = [
        {
            "source_id": 1,
            "ra": 10.0,
            "dec": 20.0,
            "parallax": 2.0,
            "parallax_error": 0.2,
            "pmra": 0.0,
            "pmra_error": 1.0,
            "pmdec": 0.0,
            "pmdec_error": 1.0,
            "classprob_dsc_combmod_star": 0.2,
        }
    ]
    result = GaiaProvider().normalize(rows)
    assert result[0].extra["secure_star"] is True


def test_legacy_query_uses_spherical_index_and_photo_z_join():
    context = TransientContext(CircleLocalization(359.9, -20.0, 10.0))
    query = LegacySurveyProvider().adql(
        context, SearchConfig(mode="full", image_recovery=False)
    )
    assert "q3c_radial_query" in query
    assert "LEFT JOIN ls_dr10.photo_z" in query
    assert "BETWEEN" not in query


def test_provider_query_uses_explicit_search_geometry():
    context = TransientContext(
        CircleLocalization(10.0, 20.0, 180.0),
        optical_ra_deg=10.1,
        optical_dec_deg=20.1,
    )
    config = SearchConfig(
        mode="full",
        search_ra_deg=10.1,
        search_dec_deg=20.1,
        search_radius_arcsec=20.0,
        image_recovery=False,
    )
    query = LegacySurveyProvider().adql(context, config)
    assert "10.10000000,20.10000000,0.00555556" in query


def test_desi_query_uses_released_coordinate_and_boolean_columns():
    context = TransientContext(CircleLocalization(10.0, 20.0, 10.0))
    query = DesiProvider().adql(
        context, SearchConfig(mode="full", image_recovery=False)
    )
    assert "mean_fiber_ra" in query
    assert "zcat_primary='t'" in query


def test_ned_current_field_names_normalize():
    candidates = NedProvider().normalize(
        [
            {
                "prefname": "Galaxy A",
                "ra": 10.0,
                "dec": 20.0,
                "z": 0.2,
                "zflag": "SPEC",
                "ptype": "G",
            }
        ]
    )
    assert candidates[0].name == "Galaxy A"
    assert candidates[0].redshift_kind == "spec"


def test_panstarrs_rejects_missing_magnitude_sentinel():
    candidates = PanStarrsProvider().normalize(
        [
            {
                "objID": 123,
                "raMean": 10.0,
                "decMean": 20.0,
                "rMeanPSFMag": -999.0,
                "rMeanKronMag": -999.0,
            }
        ]
    )
    assert candidates[0].magnitude_r is None
    assert candidates[0].morphology == "compact"
