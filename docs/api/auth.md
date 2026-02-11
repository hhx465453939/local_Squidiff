# Auth API

Lightweight intranet authentication based on local SQLite storage.

## POST `/api/auth/register`

Register a new user and create a login session.

Request JSON:
- `username` (string, 3-32, letters/numbers/`_`/`.`/`-`)
- `password` (string, min 8)

Response:
- `200 OK` -> `{ access_token, token_type: "bearer", user }`
- `409 Conflict` -> username already exists

## POST `/api/auth/login`

Login with existing credentials.

Request JSON:
- `username`
- `password`

Response:
- `200 OK` -> `{ access_token, token_type: "bearer", user }`
- `401 Unauthorized` -> invalid username or password

## GET `/api/auth/me`

Get current user by Bearer token.

Headers:
- `Authorization: Bearer <access_token>`

Response:
- `200 OK` -> `{ user }`
- `401 Unauthorized` -> invalid/expired/missing token

## POST `/api/auth/logout`

Revoke current session token.

Headers:
- `Authorization: Bearer <access_token>`

Response:
- `200 OK` -> `{ status: "ok" }`

## GET `/api/auth/user-guide`

Render frontend user guide in browser-friendly HTML.

Query:
- `raw` (bool, optional, default `false`)
  - `false`: return rendered HTML page
  - `true`: return original markdown file stream

Response:
- `200 OK` -> rendered HTML (default) or markdown stream (`raw=true`)
