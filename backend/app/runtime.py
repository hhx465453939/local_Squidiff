from __future__ import annotations

from pathlib import Path

from .core.config import settings
from .services.job_queue import JobQueue
from .services.auth_service import AuthService
from .services.squidiff_runner import SquidiffRunner
from .storage.state_manager import JsonStateStore

store = JsonStateStore(settings.state_dir)
runner = SquidiffRunner(settings.repo_root)
auth_service = AuthService(
    db_path=Path(settings.auth_db_path),
    session_ttl_hours=settings.auth_session_ttl_hours,
)
job_queue = JobQueue(store=store, runner=runner, auth_service=auth_service)
