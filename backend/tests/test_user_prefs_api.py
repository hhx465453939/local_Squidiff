from __future__ import annotations

from uuid import uuid4

from fastapi.testclient import TestClient

from backend.app.main import app


def _register_and_token(client: TestClient) -> str:
    username = f"user_{uuid4().hex[:8]}"
    password = "labflow-pass-123"
    resp = client.post(
        "/api/auth/register",
        json={"username": username, "password": password},
    )
    assert resp.status_code == 200
    return str(resp.json()["access_token"])


def test_scheduler_pref_default_is_parallel() -> None:
    client = TestClient(app)
    token = _register_and_token(client)

    resp = client.get(
        "/api/user-prefs/scheduler",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert resp.status_code == 200
    payload = resp.json()
    assert payload["mode"] == "parallel"
    assert payload["max_concurrent_jobs"] == 3


def test_scheduler_pref_can_switch_to_parallel() -> None:
    client = TestClient(app)
    token = _register_and_token(client)

    update = client.put(
        "/api/user-prefs/scheduler",
        json={"mode": "parallel"},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert update.status_code == 200
    assert update.json()["mode"] == "parallel"
    assert update.json()["max_concurrent_jobs"] == 3

    read_back = client.get(
        "/api/user-prefs/scheduler",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert read_back.status_code == 200
    assert read_back.json()["mode"] == "parallel"
