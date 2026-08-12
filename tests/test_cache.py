import json
import time

from host_copilot.cache import QueryCache


def test_query_cache_round_trip_and_release_key(tmp_path):
    cache = QueryCache(tmp_path)
    parameters = {"ra": 10.0, "dec": 20.0, "z_max": 0.1}
    first_key = cache.make_key("regalade", "v1", parameters)
    second_key = cache.make_key("regalade", "v2", parameters)
    assert first_key != second_key
    cache.put("regalade", first_key, [{"name": "host"}])
    entry = cache.get("regalade", first_key, 60, 3600)
    assert entry is not None
    assert entry.fresh
    assert entry.payload == [{"name": "host"}]


def test_query_cache_marks_old_entry_stale_but_usable(tmp_path):
    cache = QueryCache(tmp_path)
    key = cache.make_key("regalade", "v1", {"ra": 10.0})
    cache.put("regalade", key, [])
    path = tmp_path / "regalade" / f"{key}.json"
    document = json.loads(path.read_text())
    document["created_at"] = time.time() - 120
    path.write_text(json.dumps(document))
    entry = cache.get("regalade", key, 60, 3600)
    assert entry is not None
    assert not entry.fresh
    assert entry.stale_usable
