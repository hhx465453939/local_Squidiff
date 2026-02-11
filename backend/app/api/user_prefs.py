from __future__ import annotations

import sqlite3
from typing import Literal

from fastapi import APIRouter, Depends, HTTPException
from pydantic import BaseModel

from ..auth import get_current_user
from ..runtime import auth_service

router = APIRouter(prefix="/api/user-prefs", tags=["user-prefs"])


class SchedulerModePayload(BaseModel):
    mode: Literal["serial", "parallel"]


def _mode_to_limit(mode: str) -> int:
    return 3 if mode == "parallel" else 1


@router.get("/scheduler")
async def get_scheduler_mode(
    user: dict[str, object] = Depends(get_current_user),
) -> dict[str, object]:
    user_id = int(user["id"])
    try:
        mode = auth_service.get_scheduler_mode(user_id)
    except sqlite3.Error as exc:
        raise HTTPException(
            status_code=503,
            detail=f"Auth database unavailable: {exc}",
        ) from exc
    return {
        "mode": mode,
        "max_concurrent_jobs": _mode_to_limit(mode),
    }


@router.put("/scheduler")
async def update_scheduler_mode(
    payload: SchedulerModePayload,
    user: dict[str, object] = Depends(get_current_user),
) -> dict[str, object]:
    user_id = int(user["id"])
    try:
        mode = auth_service.set_scheduler_mode(user_id, payload.mode)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except sqlite3.Error as exc:
        raise HTTPException(
            status_code=503,
            detail=f"Auth database unavailable: {exc}",
        ) from exc
    return {
        "mode": mode,
        "max_concurrent_jobs": _mode_to_limit(mode),
    }
