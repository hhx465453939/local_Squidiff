from __future__ import annotations

import html
import re
import sqlite3
from pathlib import Path
from typing import Any

from fastapi import APIRouter, Depends, HTTPException, Query
from fastapi.responses import FileResponse, HTMLResponse, Response
from pydantic import BaseModel, Field

from ..auth import get_current_token, get_current_user
from ..core.config import settings
from ..runtime import auth_service

router = APIRouter(prefix="/api/auth", tags=["auth"])

USERNAME_PATTERN = re.compile(r"^[A-Za-z0-9_.-]{3,32}$")


class AuthPayload(BaseModel):
    username: str = Field(min_length=3, max_length=32)
    password: str = Field(min_length=8, max_length=128)


def _find_user_guide_path() -> Path:
    docs_dir = settings.repo_root / "docs"
    if not docs_dir.exists():
        raise HTTPException(status_code=404, detail="User guide directory not found")

    preferred_files = [
        docs_dir / "LabFlow前端用户操作说明.md",
        docs_dir / "LabFlow用户操作说明.md",
    ]
    for file_path in preferred_files:
        if file_path.exists() and file_path.is_file():
            return file_path

    for pattern in ("LabFlow*.md", "*user*guide*.md", "*.md"):
        matches = sorted(docs_dir.glob(pattern))
        if matches:
            return matches[0]

    raise HTTPException(status_code=404, detail="User guide not found")


def _read_markdown_text(path: Path) -> str:
    for encoding in ("utf-8-sig", "utf-8", "gb18030"):
        try:
            return path.read_text(encoding=encoding)
        except UnicodeDecodeError:
            continue
    raise HTTPException(status_code=500, detail="User guide encoding is not supported")


def _sanitize_href(href: str) -> str:
    normalized = href.strip()
    if normalized.startswith(("http://", "https://", "/", "#")):
        return normalized
    return "#"


def _render_inline_markdown(text: str) -> str:
    code_tokens: dict[str, str] = {}

    def _replace_code(match: re.Match[str]) -> str:
        token = f"@@CODE_{len(code_tokens)}@@"
        code_tokens[token] = f"<code>{html.escape(match.group(1), quote=False)}</code>"
        return token

    text = re.sub(r"`([^`]+)`", _replace_code, text)
    text = html.escape(text, quote=False)

    def _replace_link(match: re.Match[str]) -> str:
        label = match.group(1)
        href = html.unescape(match.group(2))
        sanitized_href = _sanitize_href(href)
        link_attrs = ""
        if sanitized_href.startswith(("http://", "https://")):
            link_attrs = ' target="_blank" rel="noreferrer"'
        safe_href = html.escape(sanitized_href, quote=True)
        return f'<a href="{safe_href}"{link_attrs}>{label}</a>'

    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", _replace_link, text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", text)
    text = re.sub(r"\*([^*]+)\*", r"<em>\1</em>", text)

    for token, code_html in code_tokens.items():
        text = text.replace(token, code_html)
    return text


def _split_table_row(line: str) -> list[str]:
    return [part.strip() for part in line.strip().strip("|").split("|")]


def _is_table_separator(line: str) -> bool:
    content = line.strip().replace(" ", "")
    if not content.startswith("|") or not content.endswith("|"):
        return False
    cells = _split_table_row(content)
    return all(
        cell.replace("-", "").replace(":", "") == "" and "-" in cell for cell in cells
    )


def _markdown_to_html(markdown_text: str) -> str:
    lines = markdown_text.splitlines()
    html_parts: list[str] = []
    paragraph_lines: list[str] = []
    list_kind: str | None = None
    in_code = False
    code_lines: list[str] = []
    i = 0

    def _flush_paragraph() -> None:
        if not paragraph_lines:
            return
        merged = " ".join(line.strip() for line in paragraph_lines).strip()
        paragraph_lines.clear()
        if merged:
            html_parts.append(f"<p>{_render_inline_markdown(merged)}</p>")

    def _flush_list() -> None:
        nonlocal list_kind
        if list_kind is not None:
            html_parts.append(f"</{list_kind}>")
            list_kind = None

    while i < len(lines):
        raw_line = lines[i]
        stripped = raw_line.strip()

        if in_code:
            if stripped.startswith("```"):
                html_parts.append(
                    "<pre><code>"
                    + html.escape("\n".join(code_lines), quote=False)
                    + "</code></pre>"
                )
                code_lines.clear()
                in_code = False
            else:
                code_lines.append(raw_line)
            i += 1
            continue

        if stripped.startswith("```"):
            _flush_paragraph()
            _flush_list()
            in_code = True
            code_lines.clear()
            i += 1
            continue

        has_table_header = (
            "|" in raw_line
            and i + 1 < len(lines)
            and lines[i + 1].strip().startswith("|")
            and _is_table_separator(lines[i + 1])
        )
        if has_table_header:
            _flush_paragraph()
            _flush_list()
            headers = _split_table_row(raw_line)
            rows: list[list[str]] = []
            i += 2
            while i < len(lines):
                row = lines[i].strip()
                if not row.startswith("|"):
                    break
                rows.append(_split_table_row(row))
                i += 1

            table_html = ["<table><thead><tr>"]
            for cell in headers:
                table_html.append(f"<th>{_render_inline_markdown(cell)}</th>")
            table_html.append("</tr></thead><tbody>")
            for row in rows:
                table_html.append("<tr>")
                for cell in row:
                    table_html.append(f"<td>{_render_inline_markdown(cell)}</td>")
                table_html.append("</tr>")
            table_html.append("</tbody></table>")
            html_parts.append("".join(table_html))
            continue

        if not stripped:
            _flush_paragraph()
            _flush_list()
            i += 1
            continue

        heading = re.match(r"^(#{1,6})\s+(.*)$", stripped)
        if heading:
            _flush_paragraph()
            _flush_list()
            level = len(heading.group(1))
            html_parts.append(
                f"<h{level}>{_render_inline_markdown(heading.group(2))}</h{level}>"
            )
            i += 1
            continue

        if stripped in {"---", "***"}:
            _flush_paragraph()
            _flush_list()
            html_parts.append("<hr />")
            i += 1
            continue

        blockquote = re.match(r"^>\s?(.*)$", stripped)
        if blockquote:
            _flush_paragraph()
            _flush_list()
            html_parts.append(
                f"<blockquote>{_render_inline_markdown(blockquote.group(1))}</blockquote>"
            )
            i += 1
            continue

        unordered = re.match(r"^[-*+]\s+(.*)$", stripped)
        ordered = re.match(r"^\d+\.\s+(.*)$", stripped)
        if unordered or ordered:
            _flush_paragraph()
            item = unordered.group(1) if unordered else ordered.group(1)
            next_list_kind = "ul" if unordered else "ol"
            if list_kind != next_list_kind:
                _flush_list()
                list_kind = next_list_kind
                html_parts.append(f"<{list_kind}>")
            html_parts.append(f"<li>{_render_inline_markdown(item)}</li>")
            i += 1
            continue

        paragraph_lines.append(raw_line)
        i += 1

    _flush_paragraph()
    _flush_list()

    if in_code and code_lines:
        html_parts.append(
            "<pre><code>"
            + html.escape("\n".join(code_lines), quote=False)
            + "</code></pre>"
        )

    return "\n".join(html_parts)


def _build_user_guide_html(markdown_text: str, title: str) -> str:
    content_html = _markdown_to_html(markdown_text)
    safe_title = html.escape(title, quote=False)
    return f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{safe_title}</title>
  <style>
    :root {{
      --bg: #f4f7fb;
      --card: #ffffff;
      --ink: #122033;
      --muted: #5a6d84;
      --line: #d6dfec;
      --accent: #0b8f9b;
      --accent-2: #2e6fe0;
      --code-bg: #0f1728;
      --code-ink: #e6edf6;
    }}
    * {{
      box-sizing: border-box;
    }}
    body {{
      margin: 0;
      color: var(--ink);
      font-family: "Segoe UI", "PingFang SC", "Microsoft YaHei", sans-serif;
      background:
        radial-gradient(circle at 10% 15%, rgba(46, 111, 224, 0.14), transparent 35%),
        radial-gradient(circle at 85% 0%, rgba(11, 143, 155, 0.17), transparent 35%),
        var(--bg);
      animation: bg-shift 14s ease-in-out infinite alternate;
    }}
    @keyframes bg-shift {{
      from {{ background-position: 0 0, 100% 0, 0 0; }}
      to {{ background-position: 6% -3%, 94% 4%, 0 0; }}
    }}
    .shell {{
      max-width: 1040px;
      margin: 0 auto;
      padding: 36px 20px 56px;
    }}
    .hero {{
      margin-bottom: 20px;
      border-radius: 18px;
      padding: 20px 24px;
      background: linear-gradient(118deg, rgba(12, 27, 61, 0.92), rgba(30, 68, 140, 0.84));
      color: #f4f9ff;
      box-shadow: 0 16px 48px rgba(8, 21, 56, 0.24);
    }}
    .hero h1 {{
      margin: 0 0 8px;
      font-size: clamp(1.3rem, 2.6vw, 2.05rem);
      letter-spacing: 0.02em;
    }}
    .hero p {{
      margin: 0;
      color: rgba(235, 246, 255, 0.88);
      font-size: 0.96rem;
    }}
    .doc {{
      border: 1px solid var(--line);
      border-radius: 18px;
      background: color-mix(in srgb, var(--card) 94%, #f4f8ff 6%);
      box-shadow: 0 14px 36px rgba(8, 31, 61, 0.09);
      padding: clamp(16px, 2.8vw, 34px);
      line-height: 1.72;
      overflow-wrap: break-word;
    }}
    h1, h2, h3, h4, h5, h6 {{
      margin: 1.3em 0 0.6em;
      line-height: 1.25;
      color: #0c3157;
    }}
    h1 {{ font-size: 1.85rem; }}
    h2 {{ font-size: 1.42rem; border-bottom: 1px solid var(--line); padding-bottom: 0.25em; }}
    h3 {{ font-size: 1.18rem; }}
    p {{
      margin: 0.56em 0;
      color: #17314c;
    }}
    ul, ol {{
      padding-left: 1.35rem;
      margin: 0.45em 0 0.95em;
    }}
    li + li {{
      margin-top: 0.24em;
    }}
    a {{
      color: var(--accent-2);
      text-decoration: none;
      border-bottom: 1px dashed rgba(46, 111, 224, 0.4);
    }}
    a:hover {{
      color: #2458b6;
      border-bottom-color: rgba(36, 88, 182, 0.72);
    }}
    hr {{
      border: 0;
      border-top: 1px solid var(--line);
      margin: 1.25rem 0;
    }}
    blockquote {{
      margin: 1.1em 0;
      border-left: 4px solid color-mix(in srgb, var(--accent) 62%, white 38%);
      background: color-mix(in srgb, #eef4ff 78%, white 22%);
      padding: 0.65em 0.95em;
      color: #2a4664;
      border-radius: 10px;
    }}
    code {{
      background: #eaf1fd;
      color: #174f8a;
      padding: 0.13em 0.42em;
      border-radius: 6px;
      font-family: "Consolas", "SFMono-Regular", Menlo, monospace;
      font-size: 0.92em;
    }}
    pre {{
      margin: 1em 0;
      padding: 12px 14px;
      border-radius: 12px;
      background: var(--code-bg);
      color: var(--code-ink);
      overflow: auto;
      box-shadow: inset 0 0 0 1px rgba(122, 158, 218, 0.2);
    }}
    pre code {{
      background: transparent;
      color: inherit;
      padding: 0;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      margin: 0.88em 0 1.25em;
      font-size: 0.95rem;
    }}
    th, td {{
      border: 1px solid var(--line);
      padding: 10px 9px;
      text-align: left;
      vertical-align: top;
    }}
    th {{
      background: #edf4ff;
      color: #183d6a;
      font-weight: 600;
    }}
    @media (max-width: 740px) {{
      .shell {{
        padding: 18px 12px 26px;
      }}
      .hero {{
        padding: 14px 15px;
      }}
      .doc {{
        border-radius: 14px;
        padding: 13px 12px 16px;
      }}
      table {{
        display: block;
        overflow-x: auto;
      }}
    }}
  </style>
</head>
<body>
  <main class="shell">
    <section class="hero">
      <h1>LabFlow User Guide</h1>
      <p>Rendered view for easier reading in browser. Use <code>?raw=true</code> for original markdown.</p>
    </section>
    <article class="doc">
      {content_html}
    </article>
  </main>
</body>
</html>
"""


def _validate_username(username: str) -> str:
    clean = username.strip()
    if not USERNAME_PATTERN.fullmatch(clean):
        raise HTTPException(
            status_code=400,
            detail=(
                "Username must be 3-32 chars and only contain letters, numbers, "
                "underscore, dot, or hyphen."
            ),
        )
    return clean


def _build_auth_response(user: dict[str, Any]) -> dict[str, Any]:
    token = auth_service.create_session(int(user["id"]))
    return {
        "access_token": token,
        "token_type": "bearer",
        "user": user,
    }


@router.post("/register")
async def register(payload: AuthPayload) -> dict[str, Any]:
    username = _validate_username(payload.username)
    try:
        user = auth_service.register(username, payload.password)
        return _build_auth_response(user)
    except sqlite3.IntegrityError:
        raise HTTPException(status_code=409, detail="Username already exists") from None
    except sqlite3.Error as exc:
        raise HTTPException(
            status_code=503,
            detail=f"Auth database unavailable: {exc}",
        ) from exc


@router.post("/login")
async def login(payload: AuthPayload) -> dict[str, Any]:
    try:
        user = auth_service.authenticate(payload.username, payload.password)
        if user is None:
            raise HTTPException(status_code=401, detail="Invalid username or password")
        return _build_auth_response(user)
    except sqlite3.Error as exc:
        raise HTTPException(
            status_code=503,
            detail=f"Auth database unavailable: {exc}",
        ) from exc


@router.post("/logout")
async def logout(token: str = Depends(get_current_token)) -> dict[str, str]:
    try:
        auth_service.revoke_session(token)
    except sqlite3.Error as exc:
        raise HTTPException(
            status_code=503,
            detail=f"Auth database unavailable: {exc}",
        ) from exc
    return {"status": "ok"}


@router.get("/me")
async def me(user: dict[str, Any] = Depends(get_current_user)) -> dict[str, Any]:
    return {"user": user}


@router.get("/user-guide", response_model=None)
async def user_guide(raw: bool = Query(default=False)) -> Response:
    guide_path = _find_user_guide_path()

    if raw:
        return FileResponse(path=guide_path, media_type="text/markdown; charset=utf-8")

    try:
        markdown_text = _read_markdown_text(guide_path)
        rendered_page = _build_user_guide_html(
            markdown_text=markdown_text, title=guide_path.stem
        )
        return HTMLResponse(content=rendered_page)
    except Exception:  # noqa: BLE001
        # Fallback to raw markdown download to avoid hard-failing guide access.
        return FileResponse(path=guide_path, media_type="text/markdown; charset=utf-8")
