# User Prefs API

用于当前登录用户的个人调度偏好（串行/并行）设置。

## GET `/api/user-prefs/scheduler`

读取当前用户的任务调度模式。

Response:
- `200 OK` ->
```json
{
  "mode": "parallel",
  "max_concurrent_jobs": 3
}
```

说明：
- `mode` 取值：`serial` 或 `parallel`
- `max_concurrent_jobs` 为该模式对应的个人并发上限
  - `serial` -> `1`
  - `parallel` -> `3`

## PUT `/api/user-prefs/scheduler`

更新当前用户的任务调度模式。

Request JSON:
- `mode` (string, required): `serial | parallel`

Response:
- `200 OK` -> 同 GET 返回结构

错误：
- `401 Unauthorized`：未登录或 token 失效
- `422 Unprocessable Entity`：`mode` 非法
