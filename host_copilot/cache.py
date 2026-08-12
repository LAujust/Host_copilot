"""Small, release-aware JSON response cache for remote catalog queries."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import time
from typing import Any


@dataclass(slots=True)
class CacheEntry:
    payload: list[dict[str, Any]]
    created_at: float
    age_seconds: float
    fresh: bool
    stale_usable: bool


class QueryCache:
    def __init__(self, directory: str | Path | None):
        self.directory = None if directory is None else Path(directory)
        if self.directory is not None:
            self.directory.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def make_key(provider: str, version: str, parameters: dict[str, Any]) -> str:
        encoded = json.dumps(
            {"provider": provider, "version": version, "parameters": parameters},
            sort_keys=True,
            separators=(",", ":"),
            default=str,
        ).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()

    def _path(self, provider: str, key: str) -> Path | None:
        if self.directory is None:
            return None
        provider_dir = self.directory / provider
        provider_dir.mkdir(parents=True, exist_ok=True)
        return provider_dir / f"{key}.json"

    def get(
        self,
        provider: str,
        key: str,
        fresh_seconds: float,
        stale_seconds: float,
    ) -> CacheEntry | None:
        path = self._path(provider, key)
        if path is None or not path.exists():
            return None
        try:
            document = json.loads(path.read_text(encoding="utf-8"))
            created_at = float(document["created_at"])
            payload = document["payload"]
            if not isinstance(payload, list):
                return None
        except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError):
            return None
        age = max(0.0, time.time() - created_at)
        return CacheEntry(
            payload=payload,
            created_at=created_at,
            age_seconds=age,
            fresh=age <= fresh_seconds,
            stale_usable=age <= stale_seconds,
        )

    def put(self, provider: str, key: str, payload: list[dict[str, Any]]) -> None:
        path = self._path(provider, key)
        if path is None:
            return
        document = {"created_at": time.time(), "payload": payload}
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(
            json.dumps(document, allow_nan=False, default=str), encoding="utf-8"
        )
        temporary.replace(path)
