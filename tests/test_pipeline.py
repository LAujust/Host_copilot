import sys
from types import SimpleNamespace

from host_copilot.models import HostCandidate, ProviderStatus, SearchConfig
from host_copilot.models import CircleLocalization, HostSearchResult, TransientContext
from host_copilot.pipeline import HostPipeline
from host_copilot.providers import ProviderResult


class FakeProvider:
    name = "regalade"
    version = "test"

    def query(self, transient, config, cache):
        candidate = HostCandidate(
            candidate_id="host",
            ra_deg=transient.localization.ra_deg,
            dec_deg=transient.localization.dec_deg,
            name="Test host",
            catalogs=[self.name],
            catalog_ids={self.name: "1"},
            redshift=0.03,
            redshift_kind="spec",
            semimajor_arcsec=5.0,
            semiminor_arcsec=3.0,
            position_angle_deg=0.0,
        )
        return ProviderResult(
            [candidate],
            ProviderStatus(
                provider=self.name,
                state="success",
                row_count=1,
                catalog_version=self.version,
            ),
        )


def test_structured_search_and_legacy_compatibility(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "host_copilot.pipeline.build_providers", lambda names: [FakeProvider()]
    )
    pipeline = HostPipeline(10.0, 20.0, 10.0, save_path=str(tmp_path))
    result = pipeline.search(
        config=SearchConfig(
            mode="quick",
            providers=("regalade",),
            image_recovery=False,
            cache_dir=tmp_path / "cache",
        )
    )
    assert result.candidates[0].name == "Test host"
    assert result.candidates[0].rank == 1
    assert result.provider_status["regalade"].state == "success"
    assert list(pipeline.galaxy_finder.reglade_df["Name"]) == ["Test host"]


def test_full_mode_surfaces_partial_provider_failure(monkeypatch, tmp_path):
    class BrokenProvider(FakeProvider):
        name = "broken"

        def query(self, transient, config, cache):
            return ProviderResult(
                [],
                ProviderStatus(
                    provider=self.name,
                    state="service_error",
                    error="offline",
                ),
            )

    monkeypatch.setattr(
        "host_copilot.pipeline.build_providers",
        lambda names: [FakeProvider(), BrokenProvider()],
    )
    pipeline = HostPipeline(10.0, 20.0, 10.0, quick=False, save_path=str(tmp_path))
    result = pipeline.search(
        config=SearchConfig(
            mode="full",
            providers=("regalade", "broken"),
            image_recovery=False,
            cache_dir=tmp_path / "cache",
        )
    )
    assert len(result.candidates) == 1
    assert not result.complete
    assert any("partial" in warning.lower() for warning in result.warnings)


def test_search_clips_candidates_to_explicit_ot_cone():
    pipeline = HostPipeline(10.0, 20.0, 180.0)
    transient = pipeline.context
    config = SearchConfig(
        mode="full",
        search_ra_deg=10.1,
        search_dec_deg=20.1,
        search_radius_arcsec=20.0,
        image_recovery=False,
    )
    inside = HostCandidate("inside", 10.1, 20.1)
    outside = HostCandidate("outside", 10.11, 20.1)

    clipped = pipeline._within_query_cone([inside, outside], transient, config)

    assert clipped == [inside]


def test_aladin_adds_ot_marker_and_error_region(monkeypatch, tmp_path):
    class FakeAladin:
        def __init__(self, **kwargs):
            self.tables = []
            self.regions = []

        def add_table(self, table, **kwargs):
            self.tables.append((table, kwargs))

        def add_graphic_overlay_from_region(self, regions):
            self.regions.extend(regions)

    monkeypatch.setitem(
        sys.modules,
        "ipyaladin",
        SimpleNamespace(Aladin=FakeAladin, EllipseError=lambda **kwargs: kwargs),
    )
    transient = TransientContext(
        CircleLocalization(10.0, 20.0, 10.0),
        optical_ra_deg=10.001,
        optical_dec_deg=20.001,
    )
    result = HostSearchResult(
        transient=transient,
        config=SearchConfig(mode="quick", image_recovery=False),
        candidates=[HostCandidate("host", 10.001, 20.001)],
    )
    pipeline = HostPipeline(localization=transient.localization, save_path=tmp_path)

    aladin = pipeline.build_aladin(result)

    assert len(aladin.tables) == 2
    assert aladin.tables[1][0]["Name"][0] == "OT"
    assert aladin.tables[1][1]["color"] == "magenta"
    assert len(aladin.regions) == 2


def test_empty_result_builds_schema_and_aladin(monkeypatch, tmp_path):
    class FakeAladin:
        def __init__(self, **kwargs):
            self.tables = []

        def add_table(self, table, **kwargs):
            self.tables.append(table)

        def add_graphic_overlay_from_region(self, regions):
            pass

    monkeypatch.setitem(
        sys.modules,
        "ipyaladin",
        SimpleNamespace(Aladin=FakeAladin, EllipseError=lambda **kwargs: kwargs),
    )
    pipeline = HostPipeline(10.0, 20.0, 10.0, save_path=tmp_path)
    result = HostSearchResult(
        transient=pipeline.context,
        config=SearchConfig(mode="quick", image_recovery=False),
    )

    table = pipeline._candidate_visualization_table(result)
    aladin = pipeline.build_aladin(result)

    assert len(table) == 0
    assert table.colnames == [
        "Name",
        "RAJ2000",
        "DEJ2000",
        "z",
        "R1",
        "R2",
        "PA",
        "sep",
        "rank",
        "relative_probability",
        "assessment",
    ]
    assert aladin.tables == []
