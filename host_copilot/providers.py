"""Remote catalog provider implementations.

Each provider returns normalized :class:`HostCandidate` instances and a
machine-readable status.  A provider failure never prevents other providers
from contributing to a full-mode search.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from io import BytesIO, StringIO
import math
import time
from typing import Any, Iterable

from astropy.io import votable
import numpy as np
import pandas as pd
import requests

from .association import finite, stable_candidate_id
from .cache import QueryCache
from .models import (
    HostCandidate,
    ProviderStatus,
    SearchConfig,
    TransientContext,
)


class ProviderSchemaError(RuntimeError):
    """Raised when a remote catalog no longer has the expected schema."""


def _bounded_request(
    method: str,
    url: str,
    timeout_seconds: float,
    **kwargs: Any,
) -> requests.Response:
    """Make at most two HTTP attempts within a shared timeout budget."""

    deadline = time.monotonic() + timeout_seconds
    last_error: Exception | None = None
    for attempt in range(2):
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            break
        attempts_left = 2 - attempt
        attempt_timeout = max(1.0, remaining / attempts_left)
        try:
            response = requests.request(method, url, timeout=attempt_timeout, **kwargs)
            response.raise_for_status()
            return response
        except requests.RequestException as exc:
            last_error = exc
    if last_error is not None:
        raise last_error
    raise requests.Timeout(f"request exceeded {timeout_seconds:.1f} seconds")


@dataclass(slots=True)
class ProviderResult:
    candidates: list[HostCandidate] = field(default_factory=list)
    status: ProviderStatus | None = None


def _scalar(value: Any) -> Any:
    if value is None or value is np.ma.masked:
        return None
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace").strip()
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if pd.isna(value):
        return None
    return value


def _records(frame: pd.DataFrame) -> list[dict[str, Any]]:
    return [
        {str(key): _scalar(value) for key, value in row.items()}
        for row in frame.to_dict(orient="records")
    ]


def _first(row: dict[str, Any], names: Iterable[str]) -> Any:
    lowered = {str(key).casefold(): value for key, value in row.items()}
    for name in names:
        value = lowered.get(name.casefold())
        if value is not None:
            return value
    return None


def _text(value: Any) -> str:
    return "" if value is None else str(value).strip()


class CatalogProvider(ABC):
    name: str
    version: str
    host_catalog: bool = True

    def parameters(
        self, transient: TransientContext, config: SearchConfig
    ) -> dict[str, Any]:
        loc = transient.localization
        return {
            "ra": round(loc.ra_deg, 8),
            "dec": round(loc.dec_deg, 8),
            "radius_arcsec": round(
                loc.enclosing_radius_arcsec + config.association_margin_arcsec, 4
            ),
            "max_rows": config.max_rows,
            "mode": config.mode,
            "z_max": config.z_max,
        }

    def query(
        self,
        transient: TransientContext,
        config: SearchConfig,
        cache: QueryCache,
    ) -> ProviderResult:
        started = time.monotonic()
        parameters = self.parameters(transient, config)
        key = cache.make_key(self.name, self.version, parameters)
        positive_ttl = config.positive_cache_days * 86400.0
        negative_ttl = config.negative_cache_days * 86400.0
        stale_ttl = config.stale_cache_days * 86400.0
        cached = cache.get(self.name, key, positive_ttl, stale_ttl)
        if cached is not None:
            fresh_ttl = negative_ttl if not cached.payload else positive_ttl
            if cached.age_seconds <= fresh_ttl:
                candidates = self.normalize(cached.payload)
                return ProviderResult(
                    candidates,
                    ProviderStatus(
                        provider=self.name,
                        state="cached",
                        elapsed_seconds=time.monotonic() - started,
                        row_count=len(cached.payload),
                        from_cache=True,
                        catalog_version=self.version,
                    ),
                )
        try:
            rows = self.fetch(transient, config)
            cache.put(self.name, key, rows)
            candidates = self.normalize(rows)
            state = "success" if rows else "empty"
            return ProviderResult(
                candidates,
                ProviderStatus(
                    provider=self.name,
                    state=state,
                    elapsed_seconds=time.monotonic() - started,
                    row_count=len(rows),
                    truncated=len(rows) >= config.max_rows,
                    catalog_version=self.version,
                    warning=(
                        f"result reached max_rows={config.max_rows}"
                        if len(rows) >= config.max_rows
                        else None
                    ),
                ),
            )
        except Exception as exc:
            if config.use_stale_on_error and cached is not None and cached.stale_usable:
                candidates = self.normalize(cached.payload)
                return ProviderResult(
                    candidates,
                    ProviderStatus(
                        provider=self.name,
                        state="cached",
                        elapsed_seconds=time.monotonic() - started,
                        row_count=len(cached.payload),
                        from_cache=True,
                        stale=True,
                        catalog_version=self.version,
                        warning=f"using stale cache after provider error: {exc}",
                    ),
                )
            if isinstance(exc, ProviderSchemaError):
                state = "schema_error"
            elif isinstance(exc, requests.Timeout):
                state = "timeout"
            else:
                state = "service_error"
            return ProviderResult(
                [],
                ProviderStatus(
                    provider=self.name,
                    state=state,
                    elapsed_seconds=time.monotonic() - started,
                    catalog_version=self.version,
                    error=f"{type(exc).__name__}: {exc}",
                ),
            )

    @abstractmethod
    def fetch(
        self, transient: TransientContext, config: SearchConfig
    ) -> list[dict[str, Any]]:
        raise NotImplementedError

    @abstractmethod
    def normalize(self, rows: list[dict[str, Any]]) -> list[HostCandidate]:
        raise NotImplementedError


class RegaladeProvider(CatalogProvider):
    name = "regalade"
    version = "J/A+A/706/A284/regalade"
    endpoint = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"

    def fetch(
        self, transient: TransientContext, config: SearchConfig
    ) -> list[dict[str, Any]]:
        radius = (
            transient.localization.enclosing_radius_arcsec
            + config.association_margin_arcsec
        )
        parameters = {
            "-source": self.version,
            "-c": (
                f"{transient.localization.ra_deg:.8f} "
                f"{transient.localization.dec_deg:+.8f}"
            ),
            "-c.rs": f"{radius:.5f}",
            "-out.max": config.max_rows,
            "-out.all": "",
            "-sort": "_r",
        }
        if config.mode == "quick":
            parameters["z"] = f"<{config.z_max}"
        response = _bounded_request(
            "GET",
            self.endpoint,
            float(config.provider_timeout_seconds),
            params=parameters,
        )
        frame = pd.read_csv(StringIO(response.text), sep="\t", comment="#")
        if frame.empty or "RAJ2000" not in frame.columns:
            return []
        frame["RAJ2000"] = pd.to_numeric(frame["RAJ2000"], errors="coerce")
        frame = frame[frame["RAJ2000"].notna()].copy()
        return _records(frame)

    def normalize(self, rows: list[dict[str, Any]]) -> list[HostCandidate]:
        normalized: list[HostCandidate] = []
        for row in rows:
            ra = finite(_first(row, ("RAJ2000", "gal_ra", "ra")))
            dec = finite(_first(row, ("DEJ2000", "gal_dec", "dec")))
            if ra is None or dec is None:
                continue
            name = _text(_first(row, ("Name", "name", "objname")))
            # ``IdCat`` is a catalog-membership bit mask, not an object ID.
            # Use the catalog name plus precise coordinates as the stable key.
            source_id = name or f"{ra:.7f},{dec:.7f}"
            redshift = finite(_first(row, ("z", "redshift")))
            redshift_error = finite(_first(row, ("e_z", "z_err")))
            z_spec = finite(_first(row, ("z_spec",)))
            redshift_kind = "spec" if z_spec == 1 else "distance"
            measurements = []
            if redshift is not None:
                measurements.append(
                    {
                        "catalog": self.name,
                        "z": redshift,
                        "z_error": redshift_error,
                        "kind": redshift_kind,
                    }
                )
            flags: list[str] = []
            if finite(_first(row, ("fRel", "f_reliability"))) == 1:
                flags.append("low_reliability")
            if finite(_first(row, ("f_simbad_zdiscrepancy",))) == 1:
                flags.append("redshift_conflict")
            normalized.append(
                HostCandidate(
                    candidate_id=stable_candidate_id(self.name, source_id, ra, dec),
                    ra_deg=ra,
                    dec_deg=dec,
                    name=name,
                    catalogs=[self.name],
                    catalog_ids={self.name: source_id},
                    semimajor_arcsec=finite(_first(row, ("R1",))),
                    semiminor_arcsec=finite(_first(row, ("R2",))),
                    position_angle_deg=finite(_first(row, ("PA",))),
                    redshift=redshift,
                    redshift_error=redshift_error,
                    redshift_kind=redshift_kind if redshift is not None else "unknown",
                    redshift_measurements=measurements,
                    magnitude_r=finite(_first(row, ("rmag", "rmagpsf"))),
                    log_stellar_mass=finite(_first(row, ("logM", "logm"))),
                    quality_flags=flags,
                    extra={
                        "distance_mpc": finite(_first(row, ("Dist", "D"))),
                        "distance_error_mpc": finite(_first(row, ("e_Dist", "D_err"))),
                    },
                )
            )
        return normalized


class TapProvider(CatalogProvider):
    tap_url = "https://datalab.noirlab.edu/tap"

    @abstractmethod
    def adql(self, transient: TransientContext, config: SearchConfig) -> str:
        raise NotImplementedError

    def fetch(
        self, transient: TransientContext, config: SearchConfig
    ) -> list[dict[str, Any]]:
        response = _bounded_request(
            "POST",
            f"{self.tap_url}/sync",
            float(config.provider_timeout_seconds),
            data={
                "REQUEST": "doQuery",
                "LANG": "ADQL",
                "FORMAT": "csv",
                "QUERY": self.adql(transient, config),
            },
        )
        if not response.text.strip():
            return []
        if response.text.lstrip().startswith("<?xml"):
            message = "TAP service returned an XML error response"
            marker = 'value="ERROR">'
            if marker in response.text:
                message = response.text.split(marker, 1)[1].split("</INFO>", 1)[0]
            raise ProviderSchemaError(message)
        return _records(pd.read_csv(StringIO(response.text)))


def _cone_clause(
    ra_column: str, dec_column: str, transient: TransientContext, config: SearchConfig
) -> str:
    loc = transient.localization
    radius_deg = (
        loc.enclosing_radius_arcsec + config.association_margin_arcsec
    ) / 3600.0
    return (
        f"'t'=q3c_radial_query({ra_column},{dec_column},"
        f"{loc.ra_deg:.8f},{loc.dec_deg:.8f},{radius_deg:.8f})"
    )


class LegacySurveyProvider(TapProvider):
    name = "legacy"
    version = "ls_dr10+photo_z"

    def adql(self, transient: TransientContext, config: SearchConfig) -> str:
        return f"""
        SELECT TOP {config.max_rows}
          t.ls_id, t.ra, t.dec, t.type, t.sersic, t.shape_r, t.shape_e1,
          t.shape_e2, t.mag_g, t.mag_r, t.mag_i, t.mag_z, t.brick_primary,
          t.parallax, t.parallax_ivar, t.pmra, t.pmdec,
          p.z_spec, p.survey, p.z_phot_median, p.z_phot_std,
          p.z_phot_l68, p.z_phot_u68
        FROM ls_dr10.tractor AS t
        LEFT JOIN ls_dr10.photo_z AS p ON t.ls_id=p.ls_id
        WHERE {_cone_clause('t.ra', 't.dec', transient, config)}
          AND t.brick_primary=1
        """

    def normalize(self, rows: list[dict[str, Any]]) -> list[HostCandidate]:
        normalized: list[HostCandidate] = []
        for row in rows:
            ra = finite(_first(row, ("ra",)))
            dec = finite(_first(row, ("dec",)))
            if ra is None or dec is None:
                continue
            source_id = _text(_first(row, ("ls_id",)))
            z_spec = finite(_first(row, ("z_spec",)))
            z_photo = finite(_first(row, ("z_phot_median",)))
            if z_spec is not None and z_spec >= 0:
                z, z_kind, z_error = z_spec, "spec", None
            elif z_photo is not None and z_photo >= 0:
                z = z_photo
                z_kind = "photo"
                z_error = finite(_first(row, ("z_phot_std",)))
            else:
                z, z_kind, z_error = None, "unknown", None
            shape_r = finite(_first(row, ("shape_r",)))
            e1 = finite(_first(row, ("shape_e1",))) or 0.0
            e2 = finite(_first(row, ("shape_e2",))) or 0.0
            ellipticity = min(0.95, math.hypot(e1, e2))
            r1 = shape_r
            r2 = (
                shape_r * (1.0 - ellipticity) / (1.0 + ellipticity)
                if shape_r is not None
                else None
            )
            pa = math.degrees(0.5 * math.atan2(e2, e1)) if shape_r else None
            morphology = _text(_first(row, ("type",))) or None
            parallax = finite(_first(row, ("parallax",)))
            parallax_ivar = finite(_first(row, ("parallax_ivar",)))
            secure_star = False
            if parallax is not None and parallax_ivar and parallax_ivar > 0:
                secure_star = abs(parallax) * math.sqrt(parallax_ivar) >= 5.0
            measurements = []
            if z is not None:
                measurements.append(
                    {
                        "catalog": self.name,
                        "z": z,
                        "z_error": z_error,
                        "kind": z_kind,
                        "survey": _text(_first(row, ("survey",))),
                    }
                )
            normalized.append(
                HostCandidate(
                    candidate_id=stable_candidate_id(self.name, source_id, ra, dec),
                    ra_deg=ra,
                    dec_deg=dec,
                    name=f"LS {source_id}" if source_id else "",
                    catalogs=[self.name],
                    catalog_ids={self.name: source_id},
                    semimajor_arcsec=r1,
                    semiminor_arcsec=r2,
                    position_angle_deg=pa,
                    redshift=z,
                    redshift_error=z_error,
                    redshift_kind=z_kind,
                    redshift_measurements=measurements,
                    magnitude_r=finite(_first(row, ("mag_r",))),
                    morphology=morphology,
                    is_star=secure_star,
                    quality_flags=["legacy_psf"] if morphology == "PSF" else [],
                    extra={
                        "z_phot_l68": finite(_first(row, ("z_phot_l68",))),
                        "z_phot_u68": finite(_first(row, ("z_phot_u68",))),
                        "secure_star": secure_star,
                    },
                )
            )
        return normalized


class DesiProvider(TapProvider):
    name = "desi"
    version = "desi_dr1.zpix"

    def adql(self, transient: TransientContext, config: SearchConfig) -> str:
        return f"""
        SELECT TOP {config.max_rows}
          targetid, mean_fiber_ra, mean_fiber_dec, z, zerr, zwarn, spectype,
          zcat_primary
        FROM desi_dr1.zpix
        WHERE {_cone_clause('mean_fiber_ra', 'mean_fiber_dec', transient, config)}
          AND zcat_primary='t' AND zwarn=0 AND spectype='GALAXY'
        """

    def normalize(self, rows: list[dict[str, Any]]) -> list[HostCandidate]:
        normalized: list[HostCandidate] = []
        for row in rows:
            ra = finite(_first(row, ("mean_fiber_ra", "target_ra", "ra")))
            dec = finite(_first(row, ("mean_fiber_dec", "target_dec", "dec")))
            z = finite(_first(row, ("z",)))
            if ra is None or dec is None or z is None:
                continue
            source_id = _text(_first(row, ("targetid",)))
            z_error = finite(_first(row, ("zerr",)))
            normalized.append(
                HostCandidate(
                    candidate_id=stable_candidate_id(self.name, source_id, ra, dec),
                    ra_deg=ra,
                    dec_deg=dec,
                    name=f"DESI {source_id}",
                    catalogs=[self.name],
                    catalog_ids={self.name: source_id},
                    redshift=z,
                    redshift_error=z_error,
                    redshift_kind="spec",
                    redshift_measurements=[
                        {
                            "catalog": self.name,
                            "z": z,
                            "z_error": z_error,
                            "kind": "spec",
                        }
                    ],
                )
            )
        return normalized


class SdssProvider(TapProvider):
    name = "sdss"
    version = "sdss_dr17.specobj"

    def adql(self, transient: TransientContext, config: SearchConfig) -> str:
        return f"""
        SELECT TOP {config.max_rows}
          specobjid, bestobjid, ra, dec, z, zerr, class, zwarning
        FROM sdss_dr17.specobj
        WHERE {_cone_clause('ra', 'dec', transient, config)}
          AND class='GALAXY' AND zwarning=0
        """

    def normalize(self, rows: list[dict[str, Any]]) -> list[HostCandidate]:
        normalized: list[HostCandidate] = []
        for row in rows:
            ra = finite(_first(row, ("ra",)))
            dec = finite(_first(row, ("dec",)))
            z = finite(_first(row, ("z",)))
            if ra is None or dec is None or z is None:
                continue
            source_id = _text(_first(row, ("specobjid", "bestobjid")))
            z_error = finite(_first(row, ("zerr",)))
            normalized.append(
                HostCandidate(
                    candidate_id=stable_candidate_id(self.name, source_id, ra, dec),
                    ra_deg=ra,
                    dec_deg=dec,
                    name=f"SDSS {source_id}",
                    catalogs=[self.name],
                    catalog_ids={self.name: source_id},
                    redshift=z,
                    redshift_error=z_error,
                    redshift_kind="spec",
                    redshift_measurements=[
                        {
                            "catalog": self.name,
                            "z": z,
                            "z_error": z_error,
                            "kind": "spec",
                        }
                    ],
                )
            )
        return normalized


class GaiaProvider(TapProvider):
    name = "gaia"
    version = "gaia_dr3.gaia_source"
    host_catalog = False

    def adql(self, transient: TransientContext, config: SearchConfig) -> str:
        return f"""
        SELECT TOP {config.max_rows}
          source_id, ra, dec, parallax, parallax_error, pmra, pmra_error,
          pmdec, pmdec_error, classprob_dsc_combmod_star
        FROM gaia_dr3.gaia_source
        WHERE {_cone_clause('ra', 'dec', transient, config)}
        """

    def normalize(self, rows: list[dict[str, Any]]) -> list[HostCandidate]:
        normalized: list[HostCandidate] = []
        for row in rows:
            ra = finite(_first(row, ("ra",)))
            dec = finite(_first(row, ("dec",)))
            if ra is None or dec is None:
                continue
            source_id = _text(_first(row, ("source_id",)))
            parallax = finite(_first(row, ("parallax",)))
            parallax_error = finite(_first(row, ("parallax_error",)))
            pmra = finite(_first(row, ("pmra",)))
            pmra_error = finite(_first(row, ("pmra_error",)))
            pmdec = finite(_first(row, ("pmdec",)))
            pmdec_error = finite(_first(row, ("pmdec_error",)))
            star_probability = finite(_first(row, ("classprob_dsc_combmod_star",)))
            astrometric_sigma = max(
                abs(parallax) / parallax_error
                if parallax is not None and parallax_error and parallax_error > 0
                else 0.0,
                abs(pmra) / pmra_error
                if pmra is not None and pmra_error and pmra_error > 0
                else 0.0,
                abs(pmdec) / pmdec_error
                if pmdec is not None and pmdec_error and pmdec_error > 0
                else 0.0,
            )
            secure_star = astrometric_sigma >= 5.0 or (
                star_probability is not None and star_probability >= 0.9
            )
            normalized.append(
                HostCandidate(
                    candidate_id=stable_candidate_id(self.name, source_id, ra, dec),
                    ra_deg=ra,
                    dec_deg=dec,
                    name=f"Gaia DR3 {source_id}",
                    catalogs=[self.name],
                    catalog_ids={self.name: source_id},
                    is_star=secure_star,
                    extra={
                        "secure_star": secure_star,
                        "astrometric_sigma": astrometric_sigma,
                        "star_probability": star_probability,
                    },
                )
            )
        return normalized


class PanStarrsProvider(CatalogProvider):
    name = "panstarrs"
    version = "panstarrs-dr2-mean"
    endpoint = "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv"

    def fetch(
        self, transient: TransientContext, config: SearchConfig
    ) -> list[dict[str, Any]]:
        loc = transient.localization
        radius = (
            loc.enclosing_radius_arcsec + config.association_margin_arcsec
        ) / 3600.0
        response = _bounded_request(
            "GET",
            self.endpoint,
            float(config.provider_timeout_seconds),
            params={
                "ra": loc.ra_deg,
                "dec": loc.dec_deg,
                "radius": radius,
                "nDetections.gt": 4,
                "pagesize": config.max_rows,
            },
        )
        if not response.text.strip():
            return []
        return _records(pd.read_csv(StringIO(response.text)))[: config.max_rows]

    def normalize(self, rows: list[dict[str, Any]]) -> list[HostCandidate]:
        normalized: list[HostCandidate] = []
        for row in rows:
            ra = finite(_first(row, ("raMean", "ramean", "ra")))
            dec = finite(_first(row, ("decMean", "decmean", "dec")))
            if ra is None or dec is None:
                continue
            source_id = _text(_first(row, ("objID", "objid", "objName")))
            name = _text(_first(row, ("objName", "objname")))
            psf = finite(_first(row, ("rMeanPSFMag", "rmeanpsfmag")))
            kron = finite(_first(row, ("rMeanKronMag", "rmeankronmag")))
            extended = psf is not None and kron is not None and psf - kron > 0.05
            normalized.append(
                HostCandidate(
                    candidate_id=stable_candidate_id(self.name, source_id, ra, dec),
                    ra_deg=ra,
                    dec_deg=dec,
                    name=name,
                    catalogs=[self.name],
                    catalog_ids={self.name: source_id},
                    magnitude_r=kron or psf,
                    morphology="extended" if extended else "compact",
                    quality_flags=[] if extended else ["panstarrs_compact"],
                )
            )
        return normalized


class NedProvider(CatalogProvider):
    name = "ned"
    version = "NED-ConeSearchByPosition-2026"
    endpoint = "https://ned.ipac.caltech.edu/NED::API/ConeSearchByPosition"

    def fetch(
        self, transient: TransientContext, config: SearchConfig
    ) -> list[dict[str, Any]]:
        loc = transient.localization
        radius_arcmin = (
            loc.enclosing_radius_arcsec + config.association_margin_arcsec
        ) / 60.0
        response = _bounded_request(
            "GET",
            self.endpoint,
            float(config.provider_timeout_seconds),
            params={
                "LON": f"{loc.ra_deg}d",
                "LAT": f"{loc.dec_deg}d",
                "CSYS": "Equatorial",
                "EQUINOX": "J2000",
                "RADIUS": radius_arcmin,
                "MAXREC": config.max_rows,
            },
        )
        table = votable.parse(BytesIO(response.content)).get_first_table().to_table()
        return _records(table.to_pandas())

    def normalize(self, rows: list[dict[str, Any]]) -> list[HostCandidate]:
        normalized: list[HostCandidate] = []
        for row in rows:
            ra = finite(_first(row, ("RA", "RA(deg)", "ra")))
            dec = finite(_first(row, ("DEC", "DEC(deg)", "dec")))
            if ra is None or dec is None:
                continue
            name = _text(
                _first(row, ("prefname", "Object Name", "Object_Name", "name"))
            )
            object_type = _text(
                _first(
                    row,
                    ("ptype", "Type", "Physical Type", "Object_Type"),
                )
            )
            if object_type.casefold() in {"star", "*", "stellar"}:
                continue
            z = finite(_first(row, ("Redshift", "z")))
            z_flag = _text(_first(row, ("zflag",)))
            z_kind = "spec" if "SPEC" in z_flag.upper() else "distance"
            source_id = name or f"{ra:.7f},{dec:.7f}"
            normalized.append(
                HostCandidate(
                    candidate_id=stable_candidate_id(self.name, source_id, ra, dec),
                    ra_deg=ra,
                    dec_deg=dec,
                    name=name,
                    catalogs=[self.name],
                    catalog_ids={self.name: source_id},
                    redshift=z,
                    redshift_kind=z_kind if z is not None else "unknown",
                    redshift_measurements=(
                        [
                            {
                                "catalog": self.name,
                                "z": z,
                                "kind": z_kind,
                                "flag": z_flag,
                            }
                        ]
                        if z is not None
                        else []
                    ),
                    morphology=object_type or None,
                )
            )
        return normalized


PROVIDERS: dict[str, type[CatalogProvider]] = {
    "regalade": RegaladeProvider,
    "legacy": LegacySurveyProvider,
    "panstarrs": PanStarrsProvider,
    "desi": DesiProvider,
    "sdss": SdssProvider,
    "ned": NedProvider,
    "gaia": GaiaProvider,
}


def build_providers(names: Iterable[str]) -> list[CatalogProvider]:
    providers: list[CatalogProvider] = []
    for name in names:
        provider_type = PROVIDERS.get(name)
        if provider_type is None:
            raise ValueError(f"unknown catalog provider: {name}")
        providers.append(provider_type())
    return providers
