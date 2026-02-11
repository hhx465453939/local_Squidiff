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
