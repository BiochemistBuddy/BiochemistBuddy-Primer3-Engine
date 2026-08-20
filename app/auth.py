import os
from dataclasses import dataclass
from typing import Annotated

import jwt
from fastapi import Depends, HTTPException
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer

ALGORITHM = "HS256"
ISSUER = "biochemistbuddy"
AUDIENCE = "biochemistbuddy-platform"
EXECUTE_SCOPE = "primer-engine:execute"

bearer_scheme = HTTPBearer(auto_error=False)


@dataclass(frozen=True)
class EngineServiceIdentity:
    service_id: str
    customer_id: str
    scopes: frozenset[str]


def _configuration() -> tuple[str, str]:
    secret = os.getenv("SERVICE_JWT_SECRET", "").strip()
    customer_id = os.getenv("DEPLOYMENT_CUSTOMER_ID", "").strip()
    if len(secret) < 32 or not customer_id:
        raise RuntimeError("Engine service authentication is not configured")
    return secret, customer_id


def decode_engine_token(token: str) -> EngineServiceIdentity:
    secret, expected_customer_id = _configuration()
    payload = jwt.decode(
        token,
        secret,
        algorithms=[ALGORITHM],
        issuer=ISSUER,
        audience=AUDIENCE,
    )
    if payload.get("token_kind") != "service":
        raise jwt.InvalidTokenError("Not a service token")
    identity = EngineServiceIdentity(
        service_id=str(payload.get("sub") or "").strip(),
        customer_id=str(payload.get("customer_id") or "").strip(),
        scopes=frozenset(str(scope) for scope in payload.get("scopes") or []),
    )
    if not identity.service_id or not identity.customer_id or not identity.scopes:
        raise jwt.InvalidTokenError("Incomplete service identity")
    if identity.customer_id != expected_customer_id:
        raise PermissionError("Service token belongs to another customer")
    if EXECUTE_SCOPE not in identity.scopes:
        raise PermissionError(f"Missing service scope: {EXECUTE_SCOPE}")
    return identity


async def authorized_engine_service(
    credentials: Annotated[HTTPAuthorizationCredentials | None, Depends(bearer_scheme)],
) -> EngineServiceIdentity:
    if credentials is None or credentials.scheme.lower() != "bearer":
        raise HTTPException(status_code=401, detail="Service authentication required")
    try:
        return decode_engine_token(credentials.credentials)
    except RuntimeError as error:
        raise HTTPException(status_code=503, detail=str(error)) from error
    except PermissionError as error:
        raise HTTPException(status_code=403, detail=str(error)) from error
    except jwt.PyJWTError as error:
        raise HTTPException(status_code=401, detail="Invalid or expired service token") from error
