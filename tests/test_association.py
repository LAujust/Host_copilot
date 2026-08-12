import pytest

from host_copilot.association import (
    apply_gaia_star_veto,
    deduplicate_candidates,
    empirical_unknown_redshift_probability,
    filter_by_redshift,
    rank_candidates,
)
from host_copilot.models import (
    CircleLocalization,
    HostCandidate,
    SearchConfig,
    TransientContext,
)


def candidate(
    candidate_id,
    ra,
    dec,
    *,
    name="galaxy",
    catalog="regalade",
    z=0.03,
    z_kind="spec",
    r1=5.0,
    r2=3.0,
):
    return HostCandidate(
        candidate_id=candidate_id,
        ra_deg=ra,
        dec_deg=dec,
        name=name,
        catalogs=[catalog],
        catalog_ids={catalog: candidate_id},
        redshift=z,
        redshift_kind=z_kind if z is not None else "unknown",
        semimajor_arcsec=r1,
        semiminor_arcsec=r2,
        position_angle_deg=0.0,
    )


def test_deduplication_merges_same_named_regalade_rows_but_not_close_pair():
    first = candidate("a", 10.0, 20.0, name="Host A", catalog="regalade")
    duplicate = candidate(
        "b", 10.0 + 1.0 / 3600.0, 20.0, name="Host A", catalog="legacy"
    )
    close_pair = candidate(
        "c", 10.0 + 2.0 / 3600.0, 20.0, name="Host B", catalog="legacy"
    )
    merged = deduplicate_candidates([first, duplicate, close_pair])
    assert len(merged) == 2
    assert set(merged[0].catalogs) == {"regalade", "legacy"}


def test_quick_requires_known_redshift_below_cut():
    config = SearchConfig(mode="quick", image_recovery=False)
    low = candidate("low", 10, 20, z=0.05)
    high = candidate("high", 10, 20, z=0.2)
    unknown = candidate("unknown", 10, 20, z=None)
    assert filter_by_redshift([low, high, unknown], config) == [low]


def test_full_retains_unknown_redshift_but_rejects_secure_high_z():
    config = SearchConfig(mode="full", z_max=0.5, image_recovery=False)
    unknown = candidate("unknown", 10, 20, z=None)
    high = candidate("high", 10, 20, z=1.0)
    assert filter_by_redshift([unknown, high], config) == [unknown]


def test_ep260321_like_candidate_is_strong_and_redshift_consistent():
    transient = TransientContext(
        CircleLocalization(149.9287, 0.4177, 10.0),
        optical_ra_deg=149.928704,
        optical_dec_deg=0.418445,
        redshift=0.0343,
    )
    host = candidate(
        "host",
        149.928696,
        0.418407,
        name="SDSS J095942.88+002506.2",
        z=0.0345446,
        r1=4.86,
        r2=3.29,
    )
    ranked = rank_candidates(
        [host], transient, SearchConfig(mode="quick", image_recovery=False)
    )
    assert ranked[0].optical_separation_arcsec == pytest.approx(0.14, abs=0.03)
    assert ranked[0].redshift_score > 0.9
    assert ranked[0].assessment == "strong_positional"
    assert ranked[0].posterior_probability is None


def test_gaia_secure_star_is_penalized_not_silently_deleted():
    host = candidate("host", 10.0, 20.0)
    star = candidate("star", 10.0, 20.0, catalog="gaia", z=None)
    star.extra["secure_star"] = True
    apply_gaia_star_veto([host], [star])
    assert host.is_star
    assert "gaia_secure_star" in host.quality_flags


def test_unknown_redshift_probability_uses_same_field_magnitude_bin():
    unknown = candidate("unknown", 10, 20, z=None)
    unknown.magnitude_r = 20.0
    population = [unknown]
    for index, redshift in enumerate((0.1, 0.2, 0.3, 0.6, 0.7, 0.8)):
        item = candidate(str(index), 10, 20, z=redshift)
        item.magnitude_r = 20.2
        population.append(item)
    probability, method = empirical_unknown_redshift_probability(
        unknown, population, 0.5
    )
    assert probability == pytest.approx(0.5)
    assert method == "field_magnitude_beta"
