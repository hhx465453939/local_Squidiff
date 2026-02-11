from __future__ import annotations

import hashlib
import hmac
import secrets
import sqlite3
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


class AuthService:
    def __init__(self, db_path: Path, session_ttl_hours: int = 12) -> None:
        self.db_path = db_path
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self.session_ttl_seconds = max(1, session_ttl_hours) * 3600
        self._ensure_schema()

    def _connect(self) -> sqlite3.Connection:
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        return conn

    def _ensure_schema(self) -> None:
        with self._connect() as conn:
            conn.executescript(
                """
                CREATE TABLE IF NOT EXISTS users (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    username TEXT NOT NULL UNIQUE,
                    username_norm TEXT NOT NULL UNIQUE,
                    password_hash TEXT NOT NULL,
                    created_at TEXT NOT NULL
                );

                CREATE TABLE IF NOT EXISTS sessions (
                    token_hash TEXT PRIMARY KEY,
                    user_id INTEGER NOT NULL,
                    created_at INTEGER NOT NULL,
                    expires_at INTEGER NOT NULL,
                    FOREIGN KEY(user_id) REFERENCES users(id) ON DELETE CASCADE
                );

                CREATE INDEX IF NOT EXISTS idx_sessions_user_id ON sessions(user_id);
                CREATE INDEX IF NOT EXISTS idx_sessions_expires_at ON sessions(expires_at);

                CREATE TABLE IF NOT EXISTS user_prefs (
                    user_id INTEGER PRIMARY KEY,
                    scheduler_mode TEXT NOT NULL DEFAULT 'parallel',
                    updated_at TEXT NOT NULL,
                    FOREIGN KEY(user_id) REFERENCES users(id) ON DELETE CASCADE
                );
                """
            )

    def _hash_password(self, password: str) -> str:
        salt = secrets.token_bytes(16)
        iterations = 200_000
        digest = hashlib.pbkdf2_hmac(
            "sha256",
            password.encode("utf-8"),
            salt,
            iterations,
        )
        return f"pbkdf2_sha256${iterations}${salt.hex()}${digest.hex()}"

    def _verify_password(self, password: str, encoded: str) -> bool:
        try:
            scheme, iter_str, salt_hex, digest_hex = encoded.split("$", 3)
            if scheme != "pbkdf2_sha256":
                return False
            iterations = int(iter_str)
            salt = bytes.fromhex(salt_hex)
            expected = bytes.fromhex(digest_hex)
        except (ValueError, TypeError):
            return False

        actual = hashlib.pbkdf2_hmac(
            "sha256",
            password.encode("utf-8"),
            salt,
            iterations,
        )
        return hmac.compare_digest(actual, expected)

    def _normalize_username(self, username: str) -> tuple[str, str]:
        clean = username.strip()
        return clean, clean.lower()

    def _row_to_user(self, row: sqlite3.Row) -> dict[str, Any]:
        return {
            "id": int(row["id"]),
            "username": str(row["username"]),
            "created_at": str(row["created_at"]),
        }

    def register(self, username: str, password: str) -> dict[str, Any]:
        clean_name, norm_name = self._normalize_username(username)
        password_hash = self._hash_password(password)
        created_at = _utc_now_iso()
        with self._connect() as conn:
            conn.execute(
                """
                INSERT INTO users (username, username_norm, password_hash, created_at)
                VALUES (?, ?, ?, ?)
                """,
                (clean_name, norm_name, password_hash, created_at),
            )
            row = conn.execute(
                "SELECT id, username, created_at FROM users WHERE username_norm = ?",
                (norm_name,),
            ).fetchone()
        if row is None:
            raise RuntimeError("User creation failed")
        return self._row_to_user(row)

    def authenticate(self, username: str, password: str) -> dict[str, Any] | None:
        _, norm_name = self._normalize_username(username)
        with self._connect() as conn:
            row = conn.execute(
                """
                SELECT id, username, password_hash, created_at
                FROM users
                WHERE username_norm = ?
                """,
                (norm_name,),
            ).fetchone()
        if row is None:
            return None
        if not self._verify_password(password, str(row["password_hash"])):
            return None
        return {
            "id": int(row["id"]),
            "username": str(row["username"]),
            "created_at": str(row["created_at"]),
        }

    def create_session(self, user_id: int) -> str:
        token = secrets.token_urlsafe(32)
        token_hash = hashlib.sha256(token.encode("utf-8")).hexdigest()
        now_ts = int(time.time())
        expires_at = now_ts + self.session_ttl_seconds
        with self._connect() as conn:
            conn.execute("DELETE FROM sessions WHERE expires_at <= ?", (now_ts,))
            conn.execute(
                """
                INSERT INTO sessions (token_hash, user_id, created_at, expires_at)
                VALUES (?, ?, ?, ?)
                """,
                (token_hash, user_id, now_ts, expires_at),
            )
        return token

    def get_user_by_token(self, token: str) -> dict[str, Any] | None:
        token_hash = hashlib.sha256(token.encode("utf-8")).hexdigest()
        now_ts = int(time.time())
        with self._connect() as conn:
            row = conn.execute(
                """
                SELECT u.id, u.username, u.created_at
                FROM sessions s
                JOIN users u ON u.id = s.user_id
                WHERE s.token_hash = ? AND s.expires_at > ?
                """,
                (token_hash, now_ts),
            ).fetchone()
        if row is None:
            return None
        return self._row_to_user(row)

    def revoke_session(self, token: str) -> None:
        token_hash = hashlib.sha256(token.encode("utf-8")).hexdigest()
        with self._connect() as conn:
            conn.execute("DELETE FROM sessions WHERE token_hash = ?", (token_hash,))

    def get_scheduler_mode(self, user_id: int) -> str:
        with self._connect() as conn:
            row = conn.execute(
                "SELECT scheduler_mode FROM user_prefs WHERE user_id = ?",
                (user_id,),
            ).fetchone()
        if row is None:
            return "parallel"
        mode = str(row["scheduler_mode"]).strip().lower()
        if mode not in {"serial", "parallel"}:
            return "parallel"
        return mode

    def set_scheduler_mode(self, user_id: int, mode: str) -> str:
        normalized = mode.strip().lower()
        if normalized not in {"serial", "parallel"}:
            raise ValueError("scheduler mode must be 'serial' or 'parallel'")
        updated_at = _utc_now_iso()
        with self._connect() as conn:
            conn.execute(
                """
                INSERT INTO user_prefs (user_id, scheduler_mode, updated_at)
                VALUES (?, ?, ?)
                ON CONFLICT(user_id) DO UPDATE SET
                    scheduler_mode=excluded.scheduler_mode,
                    updated_at=excluded.updated_at
                """,
                (user_id, normalized, updated_at),
            )
        return normalized
