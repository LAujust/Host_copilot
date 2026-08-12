from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd
import pytest

from host_copilot import HostCandidate, HostSearchResult, ProviderStatus


SCRIPT_PATH = (
    Path(__file__).resolve().parents[1] / "examples" / "process_ep_transients_modes.py"
)
SPEC = importlib.util.spec_from_file_location(
    "process_ep_transients_modes", SCRIPT_PATH
)
assert SPEC is not None and SPEC.loader is not None
BATCH = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BATCH)


@pytest.mark.parametrize(
    ("value", "expected"),
    [(10, 10.0), ("10''", 10.0), ('10"', 10.0), ("3'", 180.0), ("0.1 deg", 360.0)],
)
def test_parse_radius_arcsec(value, expected):
    assert BATCH.parse_radius_arcsec(value) == pytest.approx(expected)


def test_transient_from_ep_row_preserves_localization_and_optical_position():
    transient = BATCH.transient_from_row(
        pd.Series(
            {
                "EP Name": "EP-test",
                "RA": 149.9,
                "Dec": 0.4,
                "r_err": "3'",
                "o_RA": 149.91,
                "o_Dec": 0.41,
                "Redshift": 0.035,
                "Classification": "GRB",
                "Obs Time": "2026-01-01 00:00:00",
            }
        )
    )
    assert transient.localization.enclosing_radius_arcsec == pytest.approx(180.0)
    assert transient.localization.confidence == pytest.approx(0.9)
    assert transient.optical_ra_deg == pytest.approx(149.91)
    assert transient.redshift == pytest.approx(0.035)
    assert transient.classification == "GRB"


def test_mode_switch_builds_expected_defaults(tmp_path):
    parser = BATCH.build_parser()
    quick_args = parser.parse_args(["--mode", "quick"])
    full_args = parser.parse_args(["--mode", "full"])

    transient = BATCH.transient_from_row(
        pd.Series({"EP Name": "EP-test", "RA": 10.0, "Dec": 20.0, "r_err": "10''"})
    )
    quick = BATCH.config_from_args(quick_args, tmp_path / "quick", transient)
    full = BATCH.config_from_args(full_args, tmp_path / "full", transient)

    assert quick.z_max == pytest.approx(0.1)
    assert quick.providers == ("regalade",)
    assert quick.image_recovery is False
    assert full.z_max == pytest.approx(0.5)
    assert "desi" in full.providers
    assert "gaia" in full.providers
    assert full.image_recovery is True


def test_search_geometry_uses_ep_error_plus_ten_without_ot(tmp_path):
    args = BATCH.build_parser().parse_args(["--mode", "full"])
    transient = BATCH.transient_from_row(
        pd.Series({"EP Name": "EP-test", "RA": 10.0, "Dec": 20.0, "r_err": "3'"})
    )
    config = BATCH.config_from_args(args, tmp_path, transient)
    assert config.query_geometry(transient) == pytest.approx((10.0, 20.0, 190.0))


def test_search_geometry_uses_ot_center_and_twenty_arcseconds(tmp_path):
    args = BATCH.build_parser().parse_args(["--mode", "full"])
    transient = BATCH.transient_from_row(
        pd.Series(
            {
                "EP Name": "EP-test",
                "RA": 10.0,
                "Dec": 20.0,
                "r_err": "3'",
                "o_RA": 10.01,
                "o_Dec": 20.02,
            }
        )
    )
    config = BATCH.config_from_args(args, tmp_path, transient)
    assert config.query_geometry(transient) == pytest.approx((10.01, 20.02, 20.0))


def test_selected_rows_preserves_original_indices():
    frame = pd.DataFrame({"name": ["a", "b", "c", "d"]})
    selected = BATCH.selected_rows(frame, start=1, limit=2)
    assert selected.index.tolist() == [1, 2]


def test_prior_summary_records_supports_chunked_runs(tmp_path):
    path = tmp_path / "summary.csv"
    pd.DataFrame([{"input_row": 2, "ep_name": "EP-one", "status": "success"}]).to_csv(
        path, index=False
    )

    assert BATCH.prior_summary_records(path, preserve=False) == {}
    prior = BATCH.prior_summary_records(path, preserve=True)
    assert prior[2]["ep_name"] == "EP-one"


def test_resume_record_requires_matching_search_geometry(tmp_path):
    args = BATCH.build_parser().parse_args(["--mode", "full"])
    transient = BATCH.transient_from_row(
        pd.Series({"EP Name": "EP-test", "RA": 10.0, "Dec": 20.0, "r_err": "10''"})
    )
    config = BATCH.config_from_args(args, tmp_path, transient)
    matching = {
        "mode": "full",
        "z_max": 0.5,
        "search_ra_deg": 10.0,
        "search_dec_deg": 20.0,
        "search_radius_arcsec": 20.0,
    }
    stale = dict(matching, search_radius_arcsec=130.0)

    assert BATCH.record_matches_search(matching, config, transient)
    assert not BATCH.record_matches_search(stale, config, transient)


def test_process_transient_writes_candidates_and_provenance(tmp_path, monkeypatch):
    class FakePipeline:
        def __init__(self, **kwargs):
            pass

        def search(self, transient, config):
            candidate = HostCandidate(
                candidate_id="test-candidate",
                ra_deg=10.001,
                dec_deg=20.001,
                name="Nearby galaxy",
                catalogs=["regalade"],
                redshift=0.03,
                redshift_kind="spec",
                separation_arcsec=2.0,
                relative_probability=0.8,
                rank=1,
                assessment="strong_positional",
            )
            return HostSearchResult(
                transient=transient,
                config=config,
                candidates=[candidate],
                provider_status={
                    "regalade": ProviderStatus(
                        provider="regalade", state="success", row_count=1
                    )
                },
                metadata={"elapsed_seconds": 0.1},
            )

    monkeypatch.setattr(BATCH, "HostPipeline", FakePipeline)
    args = BATCH.build_parser().parse_args(["--mode", "quick"])
    row = pd.Series(
        {
            "EP Name": "EP-test",
            "RA": 10.0,
            "Dec": 20.0,
            "r_err": "10''",
        }
    )
    record = BATCH.process_transient(
        row, input_row=2, args=args, output_root=tmp_path, directory_name="EP-test"
    )

    assert record["status"] == "success"
    assert record["candidate_count"] == 1
    assert record["top_name"] == "Nearby galaxy"
    assert record["search_center_source"] == "EP"
    assert record["search_radius_arcsec"] == pytest.approx(20.0)
    assert (tmp_path / "EP-test" / "host_candidates.csv").exists()
    assert (tmp_path / "EP-test" / "search_result.json").exists()
    assert (tmp_path / "EP-test" / "batch_record.json").exists()
