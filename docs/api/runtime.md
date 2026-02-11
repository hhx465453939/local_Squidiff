# Runtime API

Runtime discovery and machine telemetry endpoints.

## GET `/api/runtime/conda-envs`

Discover available `conda.bat` paths and Conda environment names.

Query:
- `conda_bat` (string, optional): explicit `conda.bat` path to query.

Response:
- `200 OK`
```json
{
  "conda_bat_candidates": ["F:/software/Miniconda3/condabin/conda.bat"],
  "conda_envs": ["base", "r-4.3"]
}
```

## GET `/api/runtime/gpu-stats`

Read GPU usage from `nvidia-smi` for frontend polling during training.

Response:
- `200 OK`
```json
{
  "available": true,
  "reason": null,
  "updated_at": "2026-02-11T08:12:34.000000+00:00",
  "gpus": [
    {
      "index": 0,
      "name": "NVIDIA GeForce RTX 4090",
      "utilization_gpu": 88.0,
      "utilization_memory": 71.0,
      "memory_used_mb": 12124.0,
      "memory_total_mb": 24564.0,
      "temperature_c": 67.0,
      "power_draw_w": 327.5,
      "power_limit_w": 450.0
    }
  ]
}
```

When `nvidia-smi` is unavailable:
- `available` is `false`
- `gpus` is empty
- `reason` contains the failure detail
