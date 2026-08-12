"""Candidate deduplication, geometry, redshift handling, and ranking."""

from __future__ import annotations

from bisect import bisect_left, bisect_right
import hashlib
import math
from typing import Iterable

import astropy.units as u
from astropy.coordinates import SkyCoord, search_around_sky
from scipy.special import ndtr

from .models import HostCandidate, SearchConfig, TransientContext


REDSHIFT_PRIORITY = {"unknown": 0, "photo": 1, "distance": 2, "spec": 3}


def finite(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def stable_candidate_id(catalog: str, source_id: str, ra: float, dec: float) -> str:
    text = f"{catalog}|{source_id}|{ra:.8f}|{dec:.8f}".encode("utf-8")
    return hashlib.sha1(text).hexdigest()[:16]


def _same_candidate(left: HostCandidate, right: HostCandidate) -> bool:
    for catalog, source_id in left.catalog_ids.items():
        if source_id and right.catalog_ids.get(catalog) == source_id:
            return True
    separation = float(left.coord.separation(right.coord).arcsec)
    left_name = left.name.strip().casefold()
    right_names = {
        right.name.strip().casefold(),
        *(x.casefold() for x in right.aliases),
    }
    if left_name and left_name in right_names and separation <= 5.0:
        return True
    # Keep close galaxy pairs distinct unless their centroids agree tightly.
    return separation <= 0.75


def _merge_candidate(target: HostCandidate, source: HostCandidate) -> None:
    for catalog in source.catalogs:
        if catalog not in target.catalogs:
            target.catalogs.append(catalog)
    target.catalog_ids.update(
        {key: value for key, value in source.catalog_ids.items() if value}
    )
    for name in [source.name, *source.aliases]:
        if name and name != target.name and name not in target.aliases:
            target.aliases.append(name)
    if not target.name and source.name:
        target.name = source.name
    for attribute in (
        "semimajor_arcsec",
        "semiminor_arcsec",
        "position_angle_deg",
        "magnitude_r",
        "log_stellar_mass",
        "morphology",
    ):
        if (
            getattr(target, attribute) is None
            and getattr(source, attribute) is not None
        ):
            setattr(target, attribute, getattr(source, attribute))
    target.redshift_measurements.extend(source.redshift_measurements)
    if (
        REDSHIFT_PRIORITY[source.redshift_kind]
        > REDSHIFT_PRIORITY[target.redshift_kind]
    ):
        target.redshift = source.redshift
        target.redshift_error = source.redshift_error
        target.redshift_kind = source.redshift_kind
    elif (
        source.redshift is not None
        and target.redshift is not None
        and abs(source.redshift - target.redshift)
        > 5.0
        * max(
            source.redshift_error or 0.001 * (1.0 + source.redshift),
            target.redshift_error or 0.001 * (1.0 + target.redshift),
        )
    ):
        if "redshift_conflict" not in target.quality_flags:
            target.quality_flags.append("redshift_conflict")
    target.is_star = target.is_star or source.is_star
    for flag in source.quality_flags:
        if flag not in target.quality_flags:
            target.quality_flags.append(flag)
    target.extra.update(
        {key: value for key, value in source.extra.items() if key not in target.extra}
    )


def deduplicate_candidates(candidates: Iterable[HostCandidate]) -> list[HostCandidate]:
    """Greedily merge coordinate-identical catalog representations.

    REGALADE entries are considered first because they usually carry the most
    useful galaxy extent and distance metadata.  The tight 0.75 arcsec default
    prevents a chain of matches from collapsing real close pairs.
    """

    provider_priority = {
        "regalade": 0,
        "legacy": 1,
        "desi": 2,
        "sdss": 3,
        "ned": 4,
        "panstarrs": 5,
        "image": 6,
    }
    ordered = sorted(
        candidates,
        key=lambda candidate: min(
            (provider_priority.get(name, 99) for name in candidate.catalogs),
            default=99,
        ),
    )
    if not ordered:
        return []

    # Build the spherical neighbor index once. The former implementation
    # constructed two SkyCoord objects for every candidate pair, which made a
    # typical full-mode field with several thousand sources prohibitively
    # slow. Five arcseconds covers both the 0.75-arcsec coordinate match and
    # the broader same-name test.
    coordinates = SkyCoord(
        [candidate.ra_deg for candidate in ordered] * u.deg,
        [candidate.dec_deg for candidate in ordered] * u.deg,
        frame="icrs",
    )
    left_indices, right_indices, _, _ = search_around_sky(
        coordinates, coordinates, 5.0 * u.arcsec
    )
    prior_neighbors: list[list[int]] = [[] for _ in ordered]
    for left, right in zip(left_indices, right_indices, strict=True):
        left_index = int(left)
        right_index = int(right)
        if right_index < left_index:
            prior_neighbors[left_index].append(right_index)

    merged: list[HostCandidate] = []
    raw_to_merged: dict[int, int] = {}
    id_to_merged: dict[tuple[str, str], int] = {}
    for raw_index, candidate in enumerate(ordered):
        possible_matches = {
            id_to_merged[(catalog, source_id)]
            for catalog, source_id in candidate.catalog_ids.items()
            if source_id and (catalog, source_id) in id_to_merged
        }
        possible_matches.update(
            raw_to_merged[neighbor] for neighbor in prior_neighbors[raw_index]
        )
        match_index = next(
            (
                index
                for index in sorted(possible_matches)
                if _same_candidate(merged[index], candidate)
            ),
            None,
        )
        if match_index is None:
            merged.append(candidate)
            match_index = len(merged) - 1
        else:
            _merge_candidate(merged[match_index], candidate)
        raw_to_merged[raw_index] = match_index
        for catalog, source_id in merged[match_index].catalog_ids.items():
            if source_id:
                id_to_merged.setdefault((catalog, source_id), match_index)
    return merged


def _ellipse_radius_toward(
    candidate: HostCandidate, position: SkyCoord
) -> float | None:
    r1 = finite(candidate.semimajor_arcsec)
    r2 = finite(candidate.semiminor_arcsec)
    pa = finite(candidate.position_angle_deg)
    if r1 is None or r2 is None or pa is None or r1 <= 0 or r2 <= 0:
        return None
    direction = float(candidate.coord.position_angle(position).deg)
    theta = math.radians(direction - pa)
    return 1.0 / math.sqrt((math.cos(theta) / r1) ** 2 + (math.sin(theta) / r2) ** 2)


def add_geometry(candidate: HostCandidate, transient: TransientContext) -> None:
    localization = transient.localization
    association_position = transient.association_position
    candidate.separation_arcsec = float(
        localization.center.separation(candidate.coord).arcsec
    )
    candidate.inside_localization = localization.contains(
        candidate.ra_deg, candidate.dec_deg
    )
    candidate.optical_separation_arcsec = (
        float(association_position.separation(candidate.coord).arcsec)
        if transient.optical_ra_deg is not None
        else None
    )

    ellipse_radius = _ellipse_radius_toward(candidate, association_position)
    association_sep = float(candidate.coord.separation(association_position).arcsec)
    if ellipse_radius is not None:
        candidate.directional_light_radius = association_sep / ellipse_radius
        point_in_footprint = association_sep <= ellipse_radius
        candidate.footprint_overlap = point_in_footprint or (
            candidate.separation_arcsec
            <= localization.enclosing_radius_arcsec + ellipse_radius
        )
    else:
        candidate.directional_light_radius = None
        candidate.footprint_overlap = candidate.inside_localization

    if transient.optical_ra_deg is not None:
        sigma = transient.optical_error_arcsec
        candidate.localization_score = math.exp(-0.5 * (association_sep / sigma) ** 2)
    else:
        candidate.localization_score = localization.relative_likelihood(
            candidate.ra_deg, candidate.dec_deg
        )

    if candidate.directional_light_radius is None:
        scale = max(1.0, localization.enclosing_radius_arcsec)
        candidate.offset_score = math.exp(-0.5 * (association_sep / scale) ** 2)
    else:
        # A two-DLR scale is deliberately broad enough for off-nuclear events.
        candidate.offset_score = math.exp(
            -0.5 * (candidate.directional_light_radius / 2.0) ** 2
        )


def redshift_probability_below(candidate: HostCandidate, z_max: float) -> float:
    if candidate.redshift is None:
        # Unknown-z sources remain eligible in full mode.  A future calibrated
        # magnitude/color prior can replace this conservative neutral value.
        return 0.5
    if candidate.redshift_kind == "photo" and candidate.redshift_error:
        return float(ndtr((z_max - candidate.redshift) / candidate.redshift_error))
    return 1.0 if candidate.redshift <= z_max else 0.0


def redshift_compatibility(
    candidate: HostCandidate, transient: TransientContext, z_max: float
) -> float:
    if transient.redshift is None:
        return redshift_probability_below(candidate, z_max)
    if candidate.redshift is None:
        return 0.5
    peculiar_floor = 0.001 * (1.0 + transient.redshift)
    sigma = math.sqrt(
        (candidate.redshift_error or peculiar_floor) ** 2
        + (transient.redshift_error or peculiar_floor) ** 2
    )
    return math.exp(-0.5 * ((candidate.redshift - transient.redshift) / sigma) ** 2)


def empirical_unknown_redshift_probability(
    candidate: HostCandidate,
    population: list[HostCandidate],
    z_max: float,
) -> tuple[float, str]:
    """Estimate ``P(z < z_max)`` from same-field, similar-magnitude objects.

    A beta(1, 1) prior prevents probabilities of exactly zero or one.  When a
    field has fewer than five suitable reference galaxies the method returns a
    neutral 0.5 and labels the fallback explicitly.
    """

    if candidate.magnitude_r is None:
        return 0.5, "neutral_no_magnitude"
    references = [
        item
        for item in population
        if item is not candidate
        and item.redshift is not None
        and item.magnitude_r is not None
        and abs(item.magnitude_r - candidate.magnitude_r) <= 1.0
        and item.redshift_kind in {"spec", "photo", "distance"}
    ]
    if len(references) < 5:
        return 0.5, "neutral_insufficient_field_calibration"
    below = sum(item.redshift <= z_max for item in references)
    return (below + 1.0) / (len(references) + 2.0), "field_magnitude_beta"


def _unknown_redshift_estimator(population: list[HostCandidate], z_max: float):
    """Build an O(log N) same-field magnitude estimator for ranking."""

    references = sorted(
        (item.magnitude_r, item.redshift)
        for item in population
        if item.redshift is not None
        and item.magnitude_r is not None
        and item.redshift_kind in {"spec", "photo", "distance"}
    )
    magnitudes = [magnitude for magnitude, _ in references]
    prefix_below = [0]
    for _, redshift in references:
        prefix_below.append(prefix_below[-1] + int(redshift <= z_max))

    def estimate(candidate: HostCandidate) -> tuple[float, str]:
        if candidate.magnitude_r is None:
            return 0.5, "neutral_no_magnitude"
        lower = bisect_left(magnitudes, candidate.magnitude_r - 1.0)
        upper = bisect_right(magnitudes, candidate.magnitude_r + 1.0)
        count = upper - lower
        if count < 5:
            return 0.5, "neutral_insufficient_field_calibration"
        below = prefix_below[upper] - prefix_below[lower]
        return (below + 1.0) / (count + 2.0), "field_magnitude_beta"

    return estimate


def filter_by_redshift(
    candidates: Iterable[HostCandidate], config: SearchConfig
) -> list[HostCandidate]:
    kept: list[HostCandidate] = []
    for candidate in candidates:
        probability = redshift_probability_below(candidate, float(config.z_max))
        if config.mode == "quick":
            if candidate.redshift is None or probability <= 0.0:
                continue
        elif probability < 0.01:
            continue
        kept.append(candidate)
    return kept


def apply_gaia_star_veto(
    candidates: list[HostCandidate], gaia_sources: list[HostCandidate]
) -> None:
    if not candidates or not gaia_sources:
        return
    gaia_coords = SkyCoord(
        [source.ra_deg for source in gaia_sources] * u.deg,
        [source.dec_deg for source in gaia_sources] * u.deg,
    )
    for candidate in candidates:
        index, separation, _ = candidate.coord.match_to_catalog_sky(gaia_coords)
        if separation.arcsec > 1.0:
            continue
        star = gaia_sources[int(index)]
        if star.extra.get("secure_star", False):
            candidate.is_star = True
            if "gaia_secure_star" not in candidate.quality_flags:
                candidate.quality_flags.append("gaia_secure_star")


def rank_candidates(
    candidates: list[HostCandidate], transient: TransientContext, config: SearchConfig
) -> list[HostCandidate]:
    """Calculate auditable score components and normalized relative weights.

    The returned ``relative_probability`` values are explicitly not calibrated
    host posteriors.  ``posterior_probability`` therefore remains ``None``.
    """

    if not candidates:
        return candidates
    _, _, query_radius_arcsec = config.query_geometry(transient)
    area_arcsec2 = math.pi * query_radius_arcsec**2
    surface_density = len(candidates) / max(area_arcsec2, 1.0)

    host_priors: list[float] = []
    for candidate in candidates:
        if candidate.log_stellar_mass is not None:
            # Clamp in log space before exponentiation. Besides expressing
            # the intended [0.01, 100] prior range, this prevents malformed
            # catalog sentinel values from overflowing here.
            exponent = max(-2.0, min(2.0, candidate.log_stellar_mass - 10.0))
            raw_prior = 10.0**exponent
        elif candidate.magnitude_r is not None:
            exponent = max(-2.0, min(2.0, -0.4 * (candidate.magnitude_r - 20.0)))
            raw_prior = 10.0**exponent
        else:
            raw_prior = 0.25
        host_priors.append(min(100.0, max(0.01, raw_prior)))
    prior_normalizer = max(host_priors)
    estimate_unknown_redshift = _unknown_redshift_estimator(
        candidates, float(config.z_max)
    )

    for candidate, host_prior in zip(candidates, host_priors, strict=True):
        add_geometry(candidate, transient)
        if transient.redshift is None and candidate.redshift is None:
            probability, method = estimate_unknown_redshift(candidate)
            candidate.redshift_score = probability
            candidate.extra["unknown_z_probability"] = probability
            candidate.extra["unknown_z_probability_method"] = method
        else:
            candidate.redshift_score = redshift_compatibility(
                candidate, transient, float(config.z_max)
            )
        sep = candidate.optical_separation_arcsec
        if sep is None:
            sep = candidate.separation_arcsec or 0.0
        expected_interlopers = math.pi * max(sep, 0.5) ** 2 * surface_density
        candidate.chance_score = math.exp(-expected_interlopers)
        candidate.host_prior_score = host_prior / prior_normalizer
        candidate.quality_score = 0.05 if candidate.is_star else 1.0
        if "redshift_conflict" in candidate.quality_flags:
            candidate.quality_score *= 0.5
        if "low_reliability" in candidate.quality_flags:
            candidate.quality_score *= 0.7
        candidate.association_score = max(
            1e-300,
            candidate.localization_score
            * candidate.offset_score
            * candidate.redshift_score
            * candidate.host_prior_score
            * candidate.chance_score
            * candidate.quality_score,
        )

    normalizer = sum(candidate.association_score for candidate in candidates)
    normalizer += config.unseen_host_prior
    ordered = sorted(
        candidates,
        key=lambda candidate: (
            -candidate.association_score,
            candidate.separation_arcsec
            if candidate.separation_arcsec is not None
            else math.inf,
            candidate.name,
        ),
    )
    for rank, candidate in enumerate(ordered, start=1):
        candidate.rank = rank
        candidate.relative_probability = candidate.association_score / normalizer
        candidate.posterior_probability = None
        if candidate.is_star:
            candidate.assessment = "stellar_contaminant"
        elif candidate.footprint_overlap and candidate.redshift_score >= 0.5:
            candidate.assessment = "strong_positional"
        elif candidate.inside_localization:
            candidate.assessment = "localization_candidate"
        elif candidate.footprint_overlap:
            candidate.assessment = "boundary_overlap"
        else:
            candidate.assessment = "outside_localization"
    return ordered
