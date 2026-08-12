"""Public data models for host-galaxy searches.

The models in this module intentionally contain no network or visualization
logic.  They provide a stable boundary between catalog providers and the host
association pipeline.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
import math
from pathlib import Path
from typing import Any, Literal

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
import pandas as pd


R90_TO_SIGMA = math.sqrt(-2.0 * math.log(0.1))


def _validate_position(ra_deg: float, dec_deg: float) -> None:
    if not math.isfinite(ra_deg) or not 0.0 <= ra_deg < 360.0:
        raise ValueError(f"RA must be finite and in [0, 360): {ra_deg!r}")
    if not math.isfinite(dec_deg) or not -90.0 <= dec_deg <= 90.0:
        raise ValueError(f"Dec must be finite and in [-90, 90]: {dec_deg!r}")


@dataclass(frozen=True, slots=True)
class CircleLocalization:
    """Circular ICRS localization.

    Parameters
    ----------
    radius_arcsec
        Containment radius.  The default confidence is 90%, matching the EP
        transient input catalog used by this project.
    confidence
        Containment probability represented by ``radius_arcsec``.
    """

    ra_deg: float
    dec_deg: float
    radius_arcsec: float
    confidence: float = 0.9

    def __post_init__(self) -> None:
        _validate_position(self.ra_deg, self.dec_deg)
        if not math.isfinite(self.radius_arcsec) or self.radius_arcsec <= 0:
            raise ValueError("radius_arcsec must be positive and finite")
        if not 0.0 < self.confidence < 1.0:
            raise ValueError("confidence must lie strictly between 0 and 1")

    @property
    def center(self) -> SkyCoord:
        return SkyCoord(self.ra_deg * u.deg, self.dec_deg * u.deg, frame="icrs")

    @property
    def enclosing_radius_arcsec(self) -> float:
        return self.radius_arcsec

    @property
    def sigma_arcsec(self) -> float:
        scale = math.sqrt(-2.0 * math.log(1.0 - self.confidence))
        return self.radius_arcsec / scale

    def separation_arcsec(self, ra_deg: float, dec_deg: float) -> float:
        return float(
            self.center.separation(
                SkyCoord(ra_deg * u.deg, dec_deg * u.deg, frame="icrs")
            ).arcsec
        )

    def contains(self, ra_deg: float, dec_deg: float) -> bool:
        return self.separation_arcsec(ra_deg, dec_deg) <= self.radius_arcsec

    def relative_likelihood(self, ra_deg: float, dec_deg: float) -> float:
        q = self.separation_arcsec(ra_deg, dec_deg) / self.sigma_arcsec
        return math.exp(-0.5 * q * q)


@dataclass(frozen=True, slots=True)
class EllipseLocalization:
    """Elliptical ICRS localization with position angle east of north."""

    ra_deg: float
    dec_deg: float
    semimajor_arcsec: float
    semiminor_arcsec: float
    position_angle_deg: float
    confidence: float = 0.9

    def __post_init__(self) -> None:
        _validate_position(self.ra_deg, self.dec_deg)
        if (
            not math.isfinite(self.semimajor_arcsec)
            or not math.isfinite(self.semiminor_arcsec)
            or self.semimajor_arcsec <= 0
            or self.semiminor_arcsec <= 0
        ):
            raise ValueError("ellipse axes must be positive and finite")
        if self.semiminor_arcsec > self.semimajor_arcsec:
            raise ValueError("semiminor_arcsec cannot exceed semimajor_arcsec")
        if not math.isfinite(self.position_angle_deg):
            raise ValueError("position_angle_deg must be finite")
        if not 0.0 < self.confidence < 1.0:
            raise ValueError("confidence must lie strictly between 0 and 1")

    @property
    def center(self) -> SkyCoord:
        return SkyCoord(self.ra_deg * u.deg, self.dec_deg * u.deg, frame="icrs")

    @property
    def enclosing_radius_arcsec(self) -> float:
        return self.semimajor_arcsec

    def normalized_radius(self, ra_deg: float, dec_deg: float) -> float:
        point = SkyCoord(ra_deg * u.deg, dec_deg * u.deg, frame="icrs")
        sep = float(self.center.separation(point).arcsec)
        pa = float(self.center.position_angle(point).deg)
        theta = math.radians(pa - self.position_angle_deg)
        along_major = sep * math.cos(theta)
        along_minor = sep * math.sin(theta)
        return math.hypot(
            along_major / self.semimajor_arcsec,
            along_minor / self.semiminor_arcsec,
        )

    def separation_arcsec(self, ra_deg: float, dec_deg: float) -> float:
        return float(
            self.center.separation(
                SkyCoord(ra_deg * u.deg, dec_deg * u.deg, frame="icrs")
            ).arcsec
        )

    def contains(self, ra_deg: float, dec_deg: float) -> bool:
        return self.normalized_radius(ra_deg, dec_deg) <= 1.0

    def relative_likelihood(self, ra_deg: float, dec_deg: float) -> float:
        scale = math.sqrt(-2.0 * math.log(1.0 - self.confidence))
        q_sigma = self.normalized_radius(ra_deg, dec_deg) * scale
        return math.exp(-0.5 * q_sigma * q_sigma)


Localization = CircleLocalization | EllipseLocalization


@dataclass(frozen=True, slots=True)
class TransientContext:
    """Transient metadata used by search and association scoring."""

    localization: Localization
    name: str = ""
    optical_ra_deg: float | None = None
    optical_dec_deg: float | None = None
    optical_error_arcsec: float = 0.5
    redshift: float | None = None
    redshift_error: float | None = None
    classification: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        one_optical = (self.optical_ra_deg is None) != (self.optical_dec_deg is None)
        if one_optical:
            raise ValueError("optical RA and Dec must be supplied together")
        if self.optical_ra_deg is not None and self.optical_dec_deg is not None:
            _validate_position(self.optical_ra_deg, self.optical_dec_deg)
        if self.optical_error_arcsec <= 0:
            raise ValueError("optical_error_arcsec must be positive")
        if self.redshift is not None and self.redshift < 0:
            raise ValueError("redshift cannot be negative")
        if self.redshift_error is not None and self.redshift_error <= 0:
            raise ValueError("redshift_error must be positive")

    @property
    def association_position(self) -> SkyCoord:
        if self.optical_ra_deg is not None and self.optical_dec_deg is not None:
            return SkyCoord(
                self.optical_ra_deg * u.deg,
                self.optical_dec_deg * u.deg,
                frame="icrs",
            )
        return self.localization.center


@dataclass(slots=True)
class SearchConfig:
    """Configuration for a quick or full host search."""

    mode: Literal["quick", "full"] = "quick"
    z_max: float | None = None
    deadline_seconds: float | None = None
    provider_timeout_seconds: float | None = None
    association_margin_arcsec: float = 120.0
    max_rows: int = 5000
    cache_dir: str | Path | None = None
    positive_cache_days: float = 30.0
    negative_cache_days: float = 1.0
    stale_cache_days: float = 180.0
    use_stale_on_error: bool = True
    providers: tuple[str, ...] | None = None
    image_recovery: bool | None = None
    unseen_host_prior: float = 0.2

    def __post_init__(self) -> None:
        if self.mode not in {"quick", "full"}:
            raise ValueError("mode must be 'quick' or 'full'")
        if self.z_max is None:
            self.z_max = 0.1 if self.mode == "quick" else 0.5
        if self.z_max <= 0:
            raise ValueError("z_max must be positive")
        if self.deadline_seconds is None:
            self.deadline_seconds = 30.0 if self.mode == "quick" else 180.0
        if self.provider_timeout_seconds is None:
            self.provider_timeout_seconds = 25.0 if self.mode == "quick" else 60.0
        if self.association_margin_arcsec < 0:
            raise ValueError("association_margin_arcsec cannot be negative")
        if self.max_rows <= 0:
            raise ValueError("max_rows must be positive")
        if self.image_recovery is None:
            self.image_recovery = self.mode == "full"
        if not 0.0 <= self.unseen_host_prior < 1.0:
            raise ValueError("unseen_host_prior must lie in [0, 1)")
        if self.providers is None:
            self.providers = (
                ("regalade",)
                if self.mode == "quick"
                else (
                    "regalade",
                    "legacy",
                    "panstarrs",
                    "desi",
                    "sdss",
                    "ned",
                    "gaia",
                )
            )


ProviderState = Literal[
    "success", "empty", "timeout", "service_error", "schema_error", "cached"
]


@dataclass(slots=True)
class ProviderStatus:
    provider: str
    state: ProviderState
    elapsed_seconds: float = 0.0
    row_count: int = 0
    from_cache: bool = False
    stale: bool = False
    truncated: bool = False
    catalog_version: str | None = None
    warning: str | None = None
    error: str | None = None


@dataclass(slots=True)
class HostCandidate:
    """One normalized, possibly multi-catalog host candidate."""

    candidate_id: str
    ra_deg: float
    dec_deg: float
    name: str = ""
    catalogs: list[str] = field(default_factory=list)
    catalog_ids: dict[str, str] = field(default_factory=dict)
    aliases: list[str] = field(default_factory=list)
    semimajor_arcsec: float | None = None
    semiminor_arcsec: float | None = None
    position_angle_deg: float | None = None
    redshift: float | None = None
    redshift_error: float | None = None
    redshift_kind: Literal["spec", "photo", "distance", "unknown"] = "unknown"
    redshift_measurements: list[dict[str, Any]] = field(default_factory=list)
    magnitude_r: float | None = None
    log_stellar_mass: float | None = None
    morphology: str | None = None
    is_star: bool = False
    quality_flags: list[str] = field(default_factory=list)
    separation_arcsec: float | None = None
    optical_separation_arcsec: float | None = None
    inside_localization: bool = False
    footprint_overlap: bool = False
    directional_light_radius: float | None = None
    localization_score: float = 0.0
    offset_score: float = 0.0
    redshift_score: float = 0.0
    host_prior_score: float = 0.0
    chance_score: float = 0.0
    quality_score: float = 1.0
    association_score: float = 0.0
    relative_probability: float = 0.0
    posterior_probability: float | None = None
    rank: int | None = None
    assessment: str = "unranked"
    extra: dict[str, Any] = field(default_factory=dict)

    @property
    def coord(self) -> SkyCoord:
        return SkyCoord(self.ra_deg * u.deg, self.dec_deg * u.deg, frame="icrs")

    def to_record(self) -> dict[str, Any]:
        record = asdict(self)
        for key in (
            "catalogs",
            "catalog_ids",
            "aliases",
            "redshift_measurements",
            "quality_flags",
            "extra",
        ):
            record[key] = json.dumps(record[key], sort_keys=True, default=str)
        return record


@dataclass(slots=True)
class HostSearchResult:
    """Structured result from :meth:`HostPipeline.search`."""

    transient: TransientContext
    config: SearchConfig
    candidates: list[HostCandidate] = field(default_factory=list)
    provider_status: dict[str, ProviderStatus] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)
    aladin: Any = None

    def to_pandas(self) -> pd.DataFrame:
        if not self.candidates:
            return pd.DataFrame(
                columns=[
                    field_name for field_name in HostCandidate.__dataclass_fields__
                ]
            )
        return pd.DataFrame([candidate.to_record() for candidate in self.candidates])

    def to_table(self) -> Table:
        return Table.from_pandas(self.to_pandas())

    def to_csv(self, path: str | Path) -> None:
        self.to_pandas().to_csv(path, index=False)

    @property
    def complete(self) -> bool:
        return all(
            status.state in {"success", "empty", "cached"}
            for status in self.provider_status.values()
        ) and not any(status.truncated for status in self.provider_status.values())
