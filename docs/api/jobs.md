# Jobs API

用于训练/预测任务提交、查询、日志查看与取消。

## POST `/api/jobs/train`

提交训练任务。

Request JSON:
- `dataset_id` (string, required)
- `prepared_dataset_id` (string, optional)
- `gene_size` (int, required, > 0)
- `output_dim` (int, required, > 0)
- `use_drug_structure` (bool, optional, default `false`)
- `control_dataset_id` (string, optional)
- `batch_size` (int, optional, default `64`)
- `lr` (float, optional, default `1e-4`)

Response:
- `200 OK` -> `{ "job": JobRecord }`

调度说明（按用户）：
- 任务会记录 `owner_user_id`（当前登录用户 ID）。
- 用户调度模式为 `serial` 时：个人并发上限 `1`。
- 用户调度模式为 `parallel` 时：个人并发上限 `3`。
- 不同用户之间并发额度互不影响。

## POST `/api/jobs/predict`

提交预测任务。

Request JSON:
- `dataset_id` (string, required)
- `model_id` (string, optional)
- `model_path` (string, optional)
- `gene_size` (int, required, > 0)
- `output_dim` (int, required, > 0)
- `use_drug_structure` (bool, optional)

Response:
- `200 OK` -> `{ "job": JobRecord }`

## GET `/api/jobs/{job_id}`

查询单个任务状态。

Response:
- `200 OK` -> `{ "job": JobRecord }`
- `404 Not Found` -> `Job not found`

`job.status` 可能值：
- `queued`
- `running`
- `success`
- `failed`
- `canceled`

`JobRecord` 额外字段：
- `owner_user_id` (int | null)：任务所属用户 ID。

## GET `/api/jobs/{job_id}/log`

读取任务日志（最多返回末尾 20000 字符）。

Response:
- `200 OK` -> `{ "log": string }`

说明：
- 当任务尚未生成日志文件（例如训练刚启动）时，也会返回 `200`，`log` 为空字符串。

Live log fallback note:
- When `job.log_path` has no content yet, backend also checks
  `backend/artifacts/jobs/{job_id}/logger/log.txt` for running training logs.

## POST `/api/jobs/{job_id}/cancel`

取消任务。

行为：
- `queued`：直接标记为 `canceled`
- `running`：发送停止信号并标记 `cancel_requested`，任务进程退出后状态更新为 `canceled`
- `success/failed/canceled`：返回 400

Response:
- `200 OK` -> `{ "job": JobRecord }`
- `400 Bad Request` -> `Job is already finished: ...`
- `404 Not Found` -> `Job not found`

Stale-running behavior:
- If job status is `running` but its process PID no longer exists (e.g. backend restarted),
  backend will auto-mark it `canceled` and return updated job.

## GET `/api/jobs`

列出任务（用于前端任务中心侧边栏）。

Response:
- `200 OK` -> `{ "items": JobRecord[] }`

## DELETE `/api/jobs/{job_id}`

删除任务记录（可选同时清理 artifacts）。

Query:
- `purge_artifacts` (bool, optional, default `true`)

约束：
- `queued/running` 任务不可删除，需先取消并等待结束。

行为：
- 删除 `jobs` 记录。
- 同时删除该 `job_id` 关联的模型与结果记录（metadata）。
- `purge_artifacts=true` 时，尝试删除 `backend/artifacts/jobs/{job_id}`。

Response:
- `200 OK` -> `{ "deleted": true, "removed_models": number, "removed_results": number }`
- `400 Bad Request` -> `Job is still active...`
- `404 Not Found` -> `Job not found`

## POST `/api/jobs/flush`

一键清盘接口（用于清理卡住的测试任务）。

Query:
- `scope` (string, optional, default `active`)
  - `active`: 仅清理 `queued/running` 任务
  - `all`: 清理全部任务
- `purge_artifacts` (bool, optional, default `true`)
- `force` (bool, optional, default `true`)

行为：
- 对运行中任务先尝试取消/终止进程树。
- 对无进程的“假 running”任务直接标记并清理。
- 级联删除关联模型/结果记录与（可选）任务 artifacts。

Response:
- `200 OK` ->
```json
{
  "scope": "active",
  "deleted_jobs": 4,
  "removed_models": 2,
  "removed_results": 2,
  "skipped_running_jobs": []
}
```
