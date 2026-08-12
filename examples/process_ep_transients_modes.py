#!/usr/bin/env python3
"""Batch-search host candidates for every transient in EP_transients.csv.

Quick mode uses REGALADE with the default ``z < 0.1`` cut.  Full mode uses
all configured catalog providers, defaults to ``z < 0.5``, retains unknown-z
galaxies, and enables image-source recovery.  Results are written below
``examples/EP_data/<mode>/`` by default.

Examples
--------
Run the rapid search for every row::

    python examples/process_ep_transients_modes.py --mode quick

Run the comprehensive search with a stricter nearby-galaxy cut::

    python examples/process_ep_transients_modes.py --mode full --z-max 0.3

Validate the input without making network requests::

    python examples/process_ep_transients_modes.py --mode full --dry-run
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile
from typing import Any, Sequence

import astropy.units as u
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
# Importing the compatibility image module initializes Matplotlib.  Keep that
# cache out of read-only home directories on batch/cluster systems.
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "host_copilot_matplotlib")
)

from host_copilot import (  # noqa: E402
    CircleLocalization,
    HostPipeline,
    SearchConfig,
    TransientContext,
)


DEFAULT_INPUT = REPO_ROOT / "EP_data" / "EP_transients.csv"
DEFAULT_OUTPUT_PARENT = REPO_ROOT / "examples" / "EP_data"
RADIUS_PATTERN = re.compile(
    r"^\s*([+]?(?:\d+(?:\.\d*)?|\.\d+))\s*"
    r"(arcsec|arcmin|deg|\"|''|″|'|′|‘|’|°)\s*$",
    re.IGNORECASE,
)

SUMMARY_COLUMNS = (
    "input_row",
    "ep_name",
    "mode",
    "z_max",
    "ep_ra_deg",
    "ep_dec_deg",
    "ep_r90_arcsec",
    "optical_ra_deg",
    "optical_dec_deg",
    "search_center_source",
    "search_ra_deg",
    "search_dec_deg",
    "search_radius_arcsec",
    "transient_redshift",
    "classification",
    "status",
    "complete",
    "candidate_count",
    "known_z_count",
    "unknown_z_count",
    "inside_localization_count",
    "top_name",
    "top_ra_deg",
    "top_dec_deg",
    "top_redshift",
    "top_redshift_kind",
    "top_separation_arcsec",
    "top_relative_probability",
    "top_assessment",
    "provider_states",
    "elapsed_seconds",
    "candidate_path",
    "result_path",
    "error",
)


def finite_float(value: Any) -> float | None:
    """Return a finite float, treating blank and pandas-missing values as None."""

    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def parse_radius_arcsec(value: Any) -> float:
    """Parse an EP ``r_err`` value as a positive number of arcseconds."""

    numeric = finite_float(value)
    if numeric is not None:
        if numeric <= 0:
            raise ValueError("r_err must be positive")
        return numeric
    if not isinstance(value, str) or not value.strip():
        raise ValueError("missing r_err")
    match = RADIUS_PATTERN.fullmatch(value)
    if match is None:
        raise ValueError(f"unsupported r_err value: {value!r}")
    number = float(match.group(1))
    if number <= 0:
        raise ValueError("r_err must be positive")
    raw_unit = match.group(2).lower()
    if raw_unit in {"arcmin", "'", "′", "‘", "’"}:
        unit = u.arcmin
    elif raw_unit in {"deg", "°"}:
        unit = u.deg
    else:
        unit = u.arcsec
    return float((number * unit).to_value(u.arcsec))


def safe_name(value: str, input_row: int) -> str:
    """Create a stable filesystem-safe directory name."""

    result = re.sub(r"[^A-Za-z0-9._+-]+", "_", value.strip()).strip("._")
    return result or f"row_{input_row}"


def transient_from_row(row: pd.Series) -> TransientContext:
    """Build a validated Host Copilot context from one EP catalog row."""

    ra = finite_float(row.get("RA"))
    dec = finite_float(row.get("Dec"))
    if ra is None or dec is None:
        raise ValueError("missing RA or Dec")
    radius = parse_radius_arcsec(row.get("r_err"))
    localization = CircleLocalization(
        ra_deg=ra,
        dec_deg=dec,
        radius_arcsec=radius,
        confidence=0.90,
    )

    optical_ra = finite_float(row.get("o_RA"))
    optical_dec = finite_float(row.get("o_Dec"))
    if (optical_ra is None) != (optical_dec is None):
        raise ValueError("only one of o_RA and o_Dec is present")

    redshift = finite_float(row.get("Redshift"))
    if redshift is not None and redshift < 0:
        raise ValueError("Redshift cannot be negative")
    raw_classification = row.get("Classification")
    classification = (
        None
        if pd.isna(raw_classification) or not str(raw_classification).strip()
        else str(raw_classification).strip()
    )
    raw_name = row.get("EP Name")
    if pd.isna(raw_name) or not str(raw_name).strip():
        raise ValueError("missing EP Name")

    return TransientContext(
        name=str(raw_name).strip(),
        localization=localization,
        optical_ra_deg=optical_ra,
        optical_dec_deg=optical_dec,
        redshift=redshift,
        classification=classification,
        metadata={
            "obs_time": None
            if pd.isna(row.get("Obs Time"))
            else str(row.get("Obs Time")),
            "priority": finite_float(row.get("Priority")),
            "sx": finite_float(row.get("Sx")),
        },
    )


def relative_path(path: Path) -> str:
    try:
        return path.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def atomic_json(path: Path, value: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True, default=str, allow_nan=False),
        encoding="utf-8",
    )
    temporary.replace(path)


def atomic_summary(path: Path, records: list[dict[str, Any]]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    pd.DataFrame(records, columns=SUMMARY_COLUMNS).to_csv(temporary, index=False)
    temporary.replace(path)


def prior_summary_records(path: Path, preserve: bool) -> dict[int, dict[str, Any]]:
    """Load prior rows when resuming or processing the catalog in chunks."""

    if not preserve or not path.exists():
        return {}
    try:
        frame = pd.read_csv(path).fillna("")
        return {
            int(record["input_row"]): record
            for record in frame.to_dict(orient="records")
            if finite_float(record.get("input_row")) is not None
        }
    except (OSError, ValueError, KeyError):
        return {}


def empty_summary(input_row: int, mode: str) -> dict[str, Any]:
    record = {column: "" for column in SUMMARY_COLUMNS}
    record.update(
        {
            "input_row": input_row,
            "mode": mode,
            "status": "invalid_input",
            "complete": False,
            "candidate_count": 0,
            "known_z_count": 0,
            "unknown_z_count": 0,
            "inside_localization_count": 0,
        }
    )
    return record


def record_matches_search(
    record: dict[str, Any], config: SearchConfig, transient: TransientContext
) -> bool:
    """Return whether a resumable result used the current mode and cone."""

    expected_ra, expected_dec, expected_radius = config.query_geometry(transient)
    actual_values = (
        finite_float(record.get("search_ra_deg")),
        finite_float(record.get("search_dec_deg")),
        finite_float(record.get("search_radius_arcsec")),
    )
    if any(value is None for value in actual_values):
        return False
    actual_ra, actual_dec, actual_radius = actual_values
    assert (
        actual_ra is not None and actual_dec is not None and actual_radius is not None
    )
    return (
        record.get("mode") == config.mode
        and math.isclose(actual_ra, expected_ra, abs_tol=1e-8)
        and math.isclose(actual_dec, expected_dec, abs_tol=1e-8)
        and math.isclose(actual_radius, expected_radius, abs_tol=1e-6)
        and math.isclose(
            finite_float(record.get("z_max")) or -1.0,
            float(config.z_max),
            abs_tol=1e-12,
        )
    )


def config_from_args(
    args: argparse.Namespace,
    output_dir: Path,
    transient: TransientContext,
) -> SearchConfig:
    providers = None
    if args.providers:
        providers = tuple(
            item.strip() for item in args.providers.split(",") if item.strip()
        )
    if transient.optical_ra_deg is not None and transient.optical_dec_deg is not None:
        search_ra = transient.optical_ra_deg
        search_dec = transient.optical_dec_deg
        search_radius = args.ot_radius
    else:
        search_ra = transient.localization.ra_deg
        search_dec = transient.localization.dec_deg
        search_radius = transient.localization.enclosing_radius_arcsec + args.ep_padding
    return SearchConfig(
        mode=args.mode,
        z_max=args.z_max,
        deadline_seconds=args.deadline,
        provider_timeout_seconds=args.provider_timeout,
        association_margin_arcsec=0.0,
        search_ra_deg=search_ra,
        search_dec_deg=search_dec,
        search_radius_arcsec=search_radius,
        cache_dir=output_dir / ".host_copilot_cache",
        providers=providers,
        image_recovery=args.image_recovery,
    )


def process_transient(
    row: pd.Series,
    input_row: int,
    args: argparse.Namespace,
    output_root: Path,
    directory_name: str,
) -> dict[str, Any]:
    """Search one transient and persist its candidates and provenance."""

    summary = empty_summary(input_row, args.mode)
    try:
        transient = transient_from_row(row)
    except Exception as exc:
        summary["error"] = f"{type(exc).__name__}: {exc}"
        return summary

    loc = transient.localization
    output_dir = output_root / directory_name
    config = config_from_args(args, output_dir, transient)
    search_ra, search_dec, search_radius = config.query_geometry(transient)
    summary.update(
        {
            "ep_name": transient.name,
            "ep_ra_deg": loc.ra_deg,
            "ep_dec_deg": loc.dec_deg,
            "ep_r90_arcsec": loc.enclosing_radius_arcsec,
            "optical_ra_deg": transient.optical_ra_deg,
            "optical_dec_deg": transient.optical_dec_deg,
            "search_center_source": (
                "OT" if transient.optical_ra_deg is not None else "EP"
            ),
            "search_ra_deg": search_ra,
            "search_dec_deg": search_dec,
            "search_radius_arcsec": search_radius,
            "transient_redshift": transient.redshift,
            "classification": transient.classification or "",
        }
    )
    candidate_path = output_dir / "host_candidates.csv"
    result_path = output_dir / "search_result.json"
    batch_record_path = output_dir / "batch_record.json"

    if args.resume and batch_record_path.exists() and candidate_path.exists():
        try:
            resumed = json.loads(batch_record_path.read_text(encoding="utf-8"))
            if record_matches_search(resumed, config, transient):
                resumed["status"] = f"resumed_{resumed['status']}"
                return resumed
        except (OSError, ValueError, KeyError):
            pass

    output_dir.mkdir(parents=True, exist_ok=True)
    summary["z_max"] = config.z_max
    pipeline = HostPipeline(
        localization=loc,
        mode=args.mode,
        save_path=str(output_dir),
    )
    try:
        result = pipeline.search(transient, config)
        result.to_csv(candidate_path)
        candidates = result.candidates
        top = candidates[0] if candidates else None
        provider_states = {
            name: status.state for name, status in result.provider_status.items()
        }
        if result.complete:
            status = "success" if candidates else "no_candidates"
        else:
            status = "partial" if candidates else "partial_no_candidates"
        summary.update(
            {
                "status": status,
                "complete": result.complete,
                "candidate_count": len(candidates),
                "known_z_count": sum(x.redshift is not None for x in candidates),
                "unknown_z_count": sum(x.redshift is None for x in candidates),
                "inside_localization_count": sum(
                    x.inside_localization for x in candidates
                ),
                "top_name": top.name if top else "",
                "top_ra_deg": top.ra_deg if top else "",
                "top_dec_deg": top.dec_deg if top else "",
                "top_redshift": top.redshift if top else "",
                "top_redshift_kind": top.redshift_kind if top else "",
                "top_separation_arcsec": (
                    top.optical_separation_arcsec
                    if top and transient.optical_ra_deg is not None
                    else top.separation_arcsec
                    if top
                    else ""
                ),
                "top_relative_probability": (top.relative_probability if top else ""),
                "top_assessment": top.assessment if top else "",
                "provider_states": json.dumps(provider_states, sort_keys=True),
                "elapsed_seconds": result.metadata.get("elapsed_seconds", ""),
                "candidate_path": relative_path(candidate_path),
                "result_path": relative_path(result_path),
            }
        )
        atomic_json(
            result_path,
            {
                "input_row": input_row,
                "transient": asdict(transient),
                "config": asdict(config),
                "complete": result.complete,
                "provider_status": {
                    name: asdict(provider_status)
                    for name, provider_status in result.provider_status.items()
                },
                "warnings": result.warnings,
                "metadata": result.metadata,
            },
        )
    except Exception as exc:
        summary.update(
            {
                "status": "error",
                "error": f"{type(exc).__name__}: {exc}",
            }
        )
    atomic_json(batch_record_path, summary)
    return summary


def positive_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("value must be positive and finite")
    return number


def nonnegative_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise argparse.ArgumentTypeError("value must be non-negative and finite")
    return number


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=("quick", "full"),
        default="quick",
        help="quick: REGALADE/z<0.1; full: all providers/z<0.5 (default: quick)",
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=None,
        help="default: examples/EP_data/<mode>",
    )
    parser.add_argument("--z-max", type=positive_float, default=None)
    parser.add_argument("--deadline", type=positive_float, default=None)
    parser.add_argument("--provider-timeout", type=positive_float, default=None)
    parser.add_argument(
        "--ep-padding",
        type=nonnegative_float,
        default=10.0,
        help="padding added to r_err when no OT is available (default: 10 arcsec)",
    )
    parser.add_argument(
        "--ot-radius",
        type=positive_float,
        default=20.0,
        help="search radius centered on an available OT (default: 20 arcsec)",
    )
    parser.add_argument(
        "--providers",
        default=None,
        help="optional comma-separated provider override",
    )
    parser.add_argument(
        "--image-recovery",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="override the mode default for image-source recovery",
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--start", type=int, default=0, help="zero-based start row")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--fail-fast", action="store_true")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="validate selected input rows without querying or writing outputs",
    )
    return parser


def selected_rows(frame: pd.DataFrame, start: int, limit: int | None) -> pd.DataFrame:
    if start < 0:
        raise ValueError("--start cannot be negative")
    if limit is not None and limit <= 0:
        raise ValueError("--limit must be positive")
    stop = None if limit is None else start + limit
    return frame.iloc[start:stop]


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    input_path = args.input.expanduser().resolve()
    if not input_path.exists():
        print(f"Input file does not exist: {input_path}", file=sys.stderr)
        return 2
    try:
        transients = selected_rows(
            pd.read_csv(input_path, encoding="utf-8-sig"), args.start, args.limit
        )
    except Exception as exc:
        print(f"Could not read input: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2

    output_root = (
        args.output_root.expanduser().resolve()
        if args.output_root is not None
        else DEFAULT_OUTPUT_PARENT / args.mode
    )
    print(
        f"Selected {len(transients)} EP transients; mode={args.mode}; "
        f"output={output_root}"
    )

    if args.dry_run:
        errors = 0
        for index, row in transients.iterrows():
            input_row = int(index) + 2
            try:
                transient = transient_from_row(row)
                config = config_from_args(args, output_root, transient)
                search_ra, search_dec, search_radius = config.query_geometry(transient)
                source = "OT" if transient.optical_ra_deg is not None else "EP"
                print(
                    f"row {input_row}: {transient.name}: valid "
                    f"(center={source} {search_ra:.6f}, {search_dec:.6f}; "
                    f"radius={search_radius:g} arcsec)"
                )
            except Exception as exc:
                errors += 1
                print(f"row {input_row}: invalid: {exc}")
        print(f"Dry run complete: {len(transients) - errors} valid, {errors} invalid")
        return 1 if errors else 0

    output_root.mkdir(parents=True, exist_ok=True)
    summary_path = output_root / "summary.csv"
    records_by_row = prior_summary_records(
        summary_path,
        preserve=args.resume or args.start > 0 or args.limit is not None,
    )
    used_directories: set[str] = set()
    error_count = 0
    for ordinal, (index, row) in enumerate(transients.iterrows(), start=1):
        input_row = int(index) + 2
        raw_name = row.get("EP Name")
        name = "" if pd.isna(raw_name) else str(raw_name).strip()
        directory_name = safe_name(name, input_row)
        if directory_name in used_directories:
            directory_name = f"{directory_name}_{input_row}"
        used_directories.add(directory_name)
        print(f"[{ordinal}/{len(transients)}] {name or f'row {input_row}'}")
        record = process_transient(row, input_row, args, output_root, directory_name)
        records_by_row[input_row] = record
        atomic_summary(
            summary_path,
            [records_by_row[key] for key in sorted(records_by_row)],
        )
        print(
            f"  status={record['status']}; " f"candidates={record['candidate_count']}"
        )
        if record["error"]:
            error_count += 1
            print(f"  error={record['error']}")
            if args.fail_fast:
                break

    print(f"Summary written to {summary_path}")
    return 1 if error_count else 0


if __name__ == "__main__":
    raise SystemExit(main())
