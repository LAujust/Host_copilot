from host_copilot.models import HostCandidate, ProviderStatus, SearchConfig
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
