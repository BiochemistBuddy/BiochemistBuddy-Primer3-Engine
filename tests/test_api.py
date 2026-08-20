from datetime import UTC, datetime, timedelta

import jwt
from fastapi.testclient import TestClient

from app.main import app

SECRET = "test-service-jwt-secret-at-least-32-characters"


def _token(*, customer_id: str = "tenant-demo", scopes: tuple[str, ...]) -> str:
    now = datetime.now(UTC)
    return jwt.encode(
        {
            "sub": "primer-design",
            "customer_id": customer_id,
            "scopes": list(scopes),
            "token_kind": "service",
            "iss": "biochemistbuddy",
            "aud": "biochemistbuddy-platform",
            "iat": int(now.timestamp()),
            "exp": int((now + timedelta(minutes=5)).timestamp()),
        },
        SECRET,
        algorithm="HS256",
    )


def test_engine_requires_authentication() -> None:
    assert TestClient(app).post("/v1/design/pcr", json={}).status_code == 401
    assert TestClient(app).post("/v1/design/internal-oligo", json={}).status_code == 401


def test_engine_rejects_cross_customer_token(monkeypatch) -> None:
    monkeypatch.setenv("SERVICE_JWT_SECRET", SECRET)
    monkeypatch.setenv("DEPLOYMENT_CUSTOMER_ID", "tenant-demo")
    token = _token(customer_id="other", scopes=("primer-engine:execute",))
    response = TestClient(app).post(
        "/v1/design/pcr",
        headers={"Authorization": f"Bearer {token}"},
        json={},
    )

    assert response.status_code == 403


def test_engine_requires_execute_scope(monkeypatch) -> None:
    monkeypatch.setenv("SERVICE_JWT_SECRET", SECRET)
    monkeypatch.setenv("DEPLOYMENT_CUSTOMER_ID", "tenant-demo")
    response = TestClient(app).post(
        "/v1/design/pcr",
        headers={"Authorization": f"Bearer {_token(scopes=('primer-design:execute',))}"},
        json={},
    )

    assert response.status_code == 403
