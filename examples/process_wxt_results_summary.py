"""Quick-mode host search for the WXT/FXT results summary.

For each input row, the script searches at the optical-transient (OT) position
with an exact 20-arcsec REGALADE cone when a named optical counterpart is
present.  Otherwise it uses the FXT position with a 20-arcsec cone, falling
back to the WXT position and a 190-arcsec cone when no refined position exists.

Examples
--------
Run the full catalog serially::

    python examples/process_wxt_results_summary.py

Run four searches concurrently::

    python examples/process_wxt_results_summary.py --workers 4

Test the first ten rows::

    python examples/process_wxt_results_summary.py --limit 10
"""

from __future__ import annotations

import argparse
import multiprocessing
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from host_copilot import CircleLocalization, HostPipeline, SearchConfig  # noqa: E402
from process_ep_transients import (  # noqa: E402
    EMPTY_CANDIDATE_COLUMNS,
    finite_float,
    nearest_candidate_fields,
    safe_directory_name,
    validate_coordinates,
    write_aladin_html,
    write_aladin_notebook,
)


INPUT_CSV = REPO_ROOT / "EP_data" / "WXT_results_summary_with_gamma_info.csv"
OUTPUT_ROOT = REPO_ROOT / "examples" / "WXT_results_summary_data"
SUMMARY_CSV = OUTPUT_ROOT / "summary.csv"
Z_CUTOUT = 0.1

FXT_SEARCH_RADIUS_ARCSEC = 20.0
WXT_SEARCH_RADIUS_ARCSEC = 190.0

# SIMBAD positions sourced from TNS.  The input CSV contains optical names but
# no OT coordinate columns, so keep the resolved coordinates explicit and
# auditable.  A future unknown optical name is rejected instead of silently
# substituting the FXT or WXT position.
OPTICAL_TRANSIENT_COORDINATES = {
    "AT 2024gsa": (191.50695417, -9.71913333),
    "AT 2024ofs": (213.95885833, -16.66081667),
    "SN 2024aihh": (345.16257917, 32.59366667),
    "SN 2025kg": (55.61828750, -22.50591111),
    "SN 2025fhm": (208.39450000, -42.80464167),
    "SN 2025wkm": (36.58605000, 37.49997778),
}

RESULT_COLUMNS = [
    "input_row",
    "fxt_position_available",
    "optical_counterpart_identified",
    "optical_name",
    "optical_ra_deg",
    "optical_dec_deg",
    "ra_deg",
    "dec_deg",
    "coordinate_source",
    "search_radius_arcsec",
    "status",
    "candidate_count",
    "nearest_name",
    "nearest_ra_deg",
    "nearest_dec_deg",
    "nearest_z",
    "nearest_sep_arcsec",
    "aladin_path",
    "aladin_notebook_path",
    "cat_table_path",
    "error",
]


def usable_position(ra_value: Any, dec_value: Any) -> tuple[float, float] | None:
    """Return a valid coordinate pair, treating (0, 0) as a missing sentinel."""

    ra = finite_float(ra_value)
    dec = finite_float(dec_value)
    if ra is None or dec is None or (ra == 0.0 and dec == 0.0):
        return None
    try:
        validate_coordinates(ra, dec, "catalog")
    except ValueError:
        return None
    return ra, dec


def select_search(
    record: dict[str, Any],
) -> tuple[float, float, str, float, bool, bool, str]:
    """Choose OT/20, FXT/20, or WXT/190 arcsec in priority order."""

    fxt_position = usable_position(record.get("fxt_ra"), record.get("fxt_dec"))
    optical_name = str(record.get("SN_name", "")).strip()
    has_optical = optical_name not in {"", "-"}
    if has_optical:
        optical_position = OPTICAL_TRANSIENT_COORDINATES.get(optical_name)
        if optical_position is None:
            raise ValueError(
                f"missing resolved coordinates for optical counterpart {optical_name!r}"
            )
        return (
            optical_position[0],
            optical_position[1],
            "optical_counterpart",
            FXT_SEARCH_RADIUS_ARCSEC,
            fxt_position is not None,
            True,
            optical_name,
        )

    if fxt_position is not None:
        return (
            fxt_position[0],
            fxt_position[1],
            "fxt_ra/fxt_dec",
            FXT_SEARCH_RADIUS_ARCSEC,
            True,
            False,
            "",
        )

    wxt_position = usable_position(record.get("wxt_ra"), record.get("wxt_dec"))
    if wxt_position is None:
        raise ValueError("missing or invalid FXT and WXT coordinates")
    return (
        wxt_position[0],
        wxt_position[1],
        "wxt_ra/wxt_dec",
        WXT_SEARCH_RADIUS_ARCSEC,
        False,
        False,
        "",
    )


def relative_output_path(path: Path) -> str:
    try:
        return path.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def empty_result(record: dict[str, Any], input_row: int) -> dict[str, Any]:
    result = dict(record)
    result.update({column: "" for column in RESULT_COLUMNS})
    result.update(
        {
            "input_row": input_row,
            "status": "invalid_input",
            "candidate_count": 0,
        }
    )
    return result


def candidates_to_dataframe(candidates: list[Any]) -> pd.DataFrame:
    """Convert structured HostCandidate results to the batch CSV schema."""

    columns = EMPTY_CANDIDATE_COLUMNS + [
        "rank",
        "relative_probability",
        "assessment",
        "inside_localization",
        "footprint_overlap",
        "redshift_score",
    ]
    records = [
        {
            "Name": candidate.name,
            "RAJ2000": candidate.ra_deg,
            "DEJ2000": candidate.dec_deg,
            "z": candidate.redshift,
            "R1": candidate.semimajor_arcsec,
            "R2": candidate.semiminor_arcsec,
            "PA": candidate.position_angle_deg,
            "sep": candidate.separation_arcsec,
            "rank": candidate.rank,
            "relative_probability": candidate.relative_probability,
            "assessment": candidate.assessment,
            "inside_localization": candidate.inside_localization,
            "footprint_overlap": candidate.footprint_overlap,
            "redshift_score": candidate.redshift_score,
        }
        for candidate in candidates
    ]
    return pd.DataFrame(records, columns=columns)


def process_source(task: dict[str, Any]) -> dict[str, Any]:
    """Worker entry point for one WXT-results row."""

    record = task["record"]
    input_row = task["input_row"]
    result = empty_result(record, input_row)
    source_name = str(record.get("name", "")).strip()
    if not source_name:
        result["error"] = "missing name"
        return result

    try:
        (
            ra,
            dec,
            coordinate_source,
            search_radius,
            has_fxt,
            has_optical,
            optical_name,
        ) = select_search(record)
    except ValueError as exc:
        result["error"] = str(exc)
        return result

    result.update(
        {
            "fxt_position_available": has_fxt,
            "optical_counterpart_identified": has_optical,
            "optical_name": optical_name,
            "optical_ra_deg": ra if has_optical else "",
            "optical_dec_deg": dec if has_optical else "",
            "ra_deg": ra,
            "dec_deg": dec,
            "coordinate_source": coordinate_source,
            "search_radius_arcsec": search_radius,
        }
    )

    output_dir = Path(task["output_root"]) / task["directory_name"]
    output_dir.mkdir(parents=True, exist_ok=True)
    cat_table_path = output_dir / "cat_table.csv"
    aladin_path = output_dir / "aladin.html"
    aladin_notebook_path = output_dir / "aladin.ipynb"

    localization = CircleLocalization(ra, dec, search_radius)
    transient_redshift = finite_float(record.get("z"))
    if transient_redshift is not None and transient_redshift <= 0:
        transient_redshift = None
    pipeline = HostPipeline(
        localization=localization,
        mode="quick",
        zcutout=Z_CUTOUT,
        save_path=str(output_dir),
        transient_name=source_name,
        transient_redshift=transient_redshift,
        optical_ra=ra if has_optical else None,
        optical_dec=dec if has_optical else None,
    )
    config = SearchConfig(
        mode="quick",
        z_max=Z_CUTOUT,
        association_margin_arcsec=0.0,
        max_rows=200,
        cache_dir=pipeline.cache_dir,
    )

    try:
        search_result = pipeline.search(config=config)
        candidate_df = candidates_to_dataframe(search_result.candidates)
    except Exception as exc:
        result.update({"status": "error", "error": str(exc)})
        return result

    provider_status = search_result.provider_status.get("regalade")
    if provider_status is None or provider_status.state not in {
        "success",
        "empty",
        "cached",
    }:
        provider_error = (
            "missing REGALADE provider status"
            if provider_status is None
            else provider_status.error or provider_status.state
        )
        result.update({"status": "error", "error": provider_error})
        return result

    candidate_df.to_csv(cat_table_path, index=False)

    # The saved views display the exact cone used for the catalog query.
    write_aladin_html(
        aladin_path, source_name, ra, dec, search_radius, candidate_df
    )
    write_aladin_notebook(
        aladin_notebook_path,
        source_name,
        ra,
        dec,
        search_radius,
        candidate_df,
    )

    result.update(
        {
            "status": "success" if not candidate_df.empty else "no_candidates",
            "candidate_count": len(candidate_df),
            "aladin_path": relative_output_path(aladin_path),
            "aladin_notebook_path": relative_output_path(aladin_notebook_path),
            "cat_table_path": relative_output_path(cat_table_path),
        }
    )
    result.update(nearest_candidate_fields(candidate_df))
    return result


def make_tasks(rows: pd.DataFrame) -> list[dict[str, Any]]:
    tasks: list[dict[str, Any]] = []
    used_directories: set[str] = set()
    for zero_based_index, row in rows.iterrows():
        input_row = zero_based_index + 2
        record = row.to_dict()
        source_name = str(record.get("name", "")).strip()
        directory_name = safe_directory_name(source_name, input_row)
        if directory_name in used_directories:
            directory_name = f"{directory_name}_{input_row}"
        used_directories.add(directory_name)
        tasks.append(
            {
                "record": record,
                "input_row": input_row,
                "directory_name": directory_name,
                "output_root": str(OUTPUT_ROOT),
            }
        )
    return tasks


def run_serial(tasks: list[dict[str, Any]]) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    total = len(tasks)
    for completed, task in enumerate(tasks, start=1):
        result = process_source(task)
        results.append(result)
        print(
            f"[{completed}/{total}] {result.get('name')}: "
            f"{result['coordinate_source'] or 'no position'}, "
            f"r={result['search_radius_arcsec'] or '-'} arcsec, "
            f"status={result['status']}, candidates={result['candidate_count']}"
        )
    return results


def run_parallel(
    tasks: list[dict[str, Any]], workers: int
) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    total = len(tasks)
    context = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(max_workers=workers, mp_context=context) as executor:
        futures = {executor.submit(process_source, task): task for task in tasks}
        for completed, future in enumerate(as_completed(futures), start=1):
            task = futures[future]
            try:
                result = future.result()
            except Exception as exc:
                result = empty_result(task["record"], task["input_row"])
                result.update({"status": "error", "error": str(exc)})
            results.append(result)
            print(
                f"[{completed}/{total}] {result.get('name')}: "
                f"{result['coordinate_source'] or 'no position'}, "
                f"r={result['search_radius_arcsec'] or '-'} arcsec, "
                f"status={result['status']}, candidates={result['candidate_count']}"
            )
    return results


def write_summary(
    results: list[dict[str, Any]], input_columns: list[str]
) -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    columns = input_columns + [
        column for column in RESULT_COLUMNS if column not in input_columns
    ]
    results.sort(key=lambda result: result["input_row"])
    temporary_path = SUMMARY_CSV.with_suffix(".csv.tmp")
    pd.DataFrame(results, columns=columns).to_csv(temporary_path, index=False)
    temporary_path.replace(SUMMARY_CSV)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Quick-search REGALADE at OT or FXT positions (20 arcsec), falling "
            "back to WXT positions (190 arcsec)."
        )
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="number of concurrent searches (default: 1)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        help="process only the first N rows (useful for testing)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.workers < 1:
        raise SystemExit("--workers must be at least 1")
    if args.limit is not None and args.limit < 1:
        raise SystemExit("--limit must be at least 1")

    rows = pd.read_csv(INPUT_CSV, encoding="utf-8-sig")
    if args.limit is not None:
        rows = rows.head(args.limit)
    tasks = make_tasks(rows)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    selections = [select_search(task["record"]) for task in tasks]
    optical_count = sum(selection[5] for selection in selections)
    fxt_count = sum(
        selection[2] == "fxt_ra/fxt_dec" for selection in selections
    )
    wxt_count = len(tasks) - optical_count - fxt_count
    print(
        f"Processing {len(tasks)} rows: {optical_count} OT/20-arcsec, "
        f"{fxt_count} FXT/20-arcsec, and {wxt_count} WXT/190-arcsec searches."
    )
    results = (
        run_serial(tasks)
        if args.workers == 1
        else run_parallel(tasks, args.workers)
    )
    write_summary(results, list(rows.columns))

    errors = sum(result["status"] == "error" for result in results)
    invalid = sum(result["status"] == "invalid_input" for result in results)
    print(f"Summary written to {SUMMARY_CSV}")
    print(f"Completed with {errors} errors and {invalid} invalid rows.")
    return 1 if errors or invalid else 0


if __name__ == "__main__":
    raise SystemExit(main())
