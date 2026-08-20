from typing import Annotated

from fastapi import Depends, FastAPI, HTTPException
from pydantic import ValidationError

from app.auth import EngineServiceIdentity, authorized_engine_service
from app.engine import design_pcr_primers
from app.models import DesignEnvelope, DesignResult

app = FastAPI(title="BiochemistBuddy Primer3 Engine", version="0.1.0")


@app.get("/healthz")
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.post("/v1/design/pcr", response_model=DesignResult)
def design_pcr(
    envelope: DesignEnvelope,
    _service: Annotated[EngineServiceIdentity, Depends(authorized_engine_service)],
) -> DesignResult:
    try:
        return design_pcr_primers(envelope.request, envelope.rule_set)
    except (KeyError, TypeError, ValidationError, ValueError) as error:
        raise HTTPException(status_code=422, detail=str(error)) from error
