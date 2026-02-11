# LabFlow MVP Debug Record

## Metadata
- Module name: labflow-mvp
- Created at: 2026-02-09
- Last updated: 2026-02-11
- Related files:
- `backend/app/main.py`
- `backend/app/api/datasets.py`
- `backend/app/api/jobs.py`
- `backend/app/api/results.py`
- `backend/app/storage/state_manager.py`
- `backend/app/services/job_queue.py`
- `backend/app/services/squidiff_runner.py`
- `backend/app/services/seurat_converter.py`
- `backend/app/services/seurat_inspector.py`
- `backend/app/services/dataset_preprocessor.py`
- `frontend/src/App.tsx`
- `frontend/src/services/api.ts`
- `backend/app/api/seurat.py`
- `backend/tests/test_seurat_api.py`
- `backend/tests/test_dataset_preprocessor.py`
- `backend/tests/test_jobs_api.py`
- `docs/api/seurat.md`
- `docs/seurat转换指南.md`
- `docs/实验室10分钟上手.md`
- `docs/UAT_Seurat_V2_检查清单.md`
- `scripts/uat_phase4_seurat_v2.py`
- `infra/docker-compose.yml`
- Dependency modules:
- `train_squidiff.py`
- `sample_squidiff.py`

## Runtime Context and Test Rules
- Runtime environment: **Windows 本机**（项目路径 `E:\Development\local_Squidiff`）；也可在 WSL2 下用 `/mnt/e/Development/local_Squidiff`。
- SSH mode (if remote): Not used.
- Remote project path (if remote): N/A
- Validation/Checkfix execution mode: 在本地 PowerShell/CMD 直接执行；Windows 下 R 须用 `cmd_conda` + `r-4.3`。
- R execution constraint (confirmed by user): R conda env must be activated via `cmd` (not PowerShell). 推荐环境：**r-4.3**（`F:\software\Miniconda3\envs\r-4.3`），包齐全、稳定。
- R config strategy: support both `.env` defaults (`LABFLOW_R_*`) and per-request frontend overrides (`r_exec_mode`, `r_conda_env`, `r_conda_bat`, `rscript_bin`).
- 示例数据（data/）：**TC.rds** = 大鼠皮下筋膜针灸/痢疾 telocytes；**coTC.rds** = 大鼠结肠针灸/痢疾 telocytes；细胞量较大，用于 500×500 流程测试。
- Windows 下 500×500 测试脚本：`scripts/run_500x500_test_windows.py`。先启动后端，再在另一终端执行 `python scripts/run_500x500_test_windows.py`；可选环境变量 `LABFLOW_BASE_URL`、`LABFLOW_R_CONDA_ENV`、`LABFLOW_R_CONDA_BAT`。
- 全流程（转换→训练→预测→可视化报告）：`scripts/run_full_train_predict_viz_windows.py`。建议后端设置 `LABFLOW_DRY_RUN=true` 以快速跑通；报告输出到 `scripts/output/full_flow_report/`（summary.json + pca_scatter.png、heatmap_top_var_genes.png）。
- **端口与 R 转换**：若 8000 被占用，可在其它端口启动后端并设 `LABFLOW_BASE_URL`。若 validate 报「conda.bat 不是内部或外部命令」，说明当前后端进程是旧代码，需**重启后端**以加载 R 转换的临时 .bat 修复（见 2026-02-11 Windows 全自动 500×500 测试条目）。
- **前端校验报 Rscript was not found**：后端进程 PATH 中无 Rscript 时（如 R 仅在 Conda 环境 r-4.3 中），校验会 400。解决：在页面上将「R 执行方式」改为 **cmd_conda**，填写 conda.bat 完整路径（如 F:\\software\\Miniconda3\\condabin\\conda.bat）和 R 环境名（如 r-4.3）。已做：后端错误文案增加 cmd_conda 说明；前端校验区增加提示、并对 400 的 detail 解析后追加操作建议。
- **前端用户说明书**：已新增 `docs/LabFlow前端用户操作说明.md`，按页面 1) 上传 2) 校验 3) Seurat 解析 4) 500x500 预处理 5) 提交训练 6) 任务轮询 7) 结果页 逐项说明每个选项与参数如何填写（含 Windows Conda R、direct/cmd_conda、常见问题）。README 第 5 节与第 9 节文档导航已引用该文档。
- **真实训练 vs dry_run**：后端设置 `LABFLOW_DRY_RUN=true` 时，训练不执行（只写占位 `model.pt`），预测用随机矩阵，图会正常生成。要得到真实训练出的模型，需**不设或关闭** `LABFLOW_DRY_RUN` 后重启后端再跑全流程；真实模型在 `backend/artifacts/jobs/<train_job_id>/checkpoints/` 下。
- **训练失败 ModuleNotFoundError: rdkit**：当 `use_drug_structure=False` 时，Squidiff 不需 rdkit。已在 `Squidiff/scrna_datasets.py` 中将 rdkit 改为在 `Drug_dose_encoder` 内按需导入，避免无药物结构时因缺 rdkit 导致训练启动失败。训练失败时后端会在 job 的 error_msg 及 train.log 中保留子进程 stderr 末尾，便于排查。
- **训练 use_drug_structure 被误传为 True**：runner 原先传 `--use_drug_structure str(False)` 即 `"False"`，argparse 解析为 True，导致脚本去读空的 control_data_path 报 OSError。已改为仅当 `params["use_drug_structure"]` 为真时才追加 `--use_drug_structure True` 与 `--control_data_path`，否则不传，使用 train_squidiff 默认 False。失败时若子进程无 stderr/stdout，则从已写入的 train.log 读末尾作为错误详情。
- **训练轮询超时但显卡仍在跑**：脚本原为固定 3600s 超时，训练超过 1 小时即报错，而 GPU 仍在训练。已在全流程脚本中增加「超时后多侧面判断」：(1) nvidia-smi 查 GPU 利用率，高于阈值则延长；(2) nvidia-smi 进程列表（--query-compute-apps 或 -q 解析）中若存在名称含 python 的进程，也视为训练可能仍在跑并延长。任一满足即延长等待（每次 30 分钟，总上限 4 小时）。环境变量：LABFLOW_TRAIN_GPU_BUSY_THRESHOLD、LABFLOW_TRAIN_EXTEND_SEC、LABFLOW_TRAIN_MAX_TOTAL_SEC。
- **训练轮询按本任务 PID 判断（精确到进程）**：不再仅看「是否有 python 在 GPU」，改为优先看**本任务训练进程**是否仍在 GPU。后端在启动训练子进程时用 Popen 获得 PID，通过 on_start(pid) 回调写入 job 的 train_pid；GET /api/jobs/{job_id} 返回该字段。脚本增加 get_gpu_pids()（nvidia-smi --query-compute-apps=pid 或 -q 解析 Process ID），超时后若 job 含 train_pid，则**仅当 train_pid in get_gpu_pids()** 时才延长；否则退化为「利用率或任意 Python 在 GPU」。这样其它程序占用 GPU 不会误触发延长。

## Context Network
- File layout
- New MVP modules added under `backend/`, `frontend/`, and `infra/`.
- Existing training/predict scripts are kept unchanged and invoked via service wrapper.
- Function call chain
- API (`datasets/jobs/results`) -> runtime singleton (`store`, `job_queue`) -> service layer (`validator/converter/runner`) -> JSON state + scripts.
- Variable/data dependencies
- Dataset upload writes raw paths into `datasets.json`.
- Job payload references dataset IDs and optional model IDs.
- Worker resolves paths then writes model/result entries into `models.json`/`results.json`.
- Data flow
- User upload -> validate/convert -> queue job -> run training/predict -> generate result assets -> query job/result APIs.
- Frontend full flow
- Upload file -> validate (with R runtime config) -> submit train job -> poll job status -> display train outputs/result assets/logs.

## Debug History
### [2026-02-09 23:xx] Bootstrap no-SQL MVP
- Problem
- Need to start implementation from PRD with no SQL requirement and minimal intranet scope.
- Root cause
- Repository only had research scripts, no service architecture.
- Solution
- Added minimal backend API, JSON state store, background worker queue, script wrappers, and frontend shell.
- Code changes (files/functions)
- Added backend runtime and API endpoints, plus file-based persistence and Docker compose deployment.
- Verification results
- `ruff check backend`: passed.
- `ruff format --check backend`: passed.
- `python -m compileall backend/app`: passed.
- `python -c "from fastapi.testclient import TestClient; ... /api/health"`: passed (`200 {"status":"ok"}`).
- `uv run pytest -q backend/tests`: blocked by dependency resolution conflict in current environment (`scanpy` vs broad `requires-python >=3.8` resolution path). Not a code syntax failure.
- `python -m pytest -q backend/tests`: blocked because local Python environment does not include `pytest`.
- Impact assessment
- Existing model scripts untouched; new code isolated under `backend/`, `frontend/`, `infra/`.

### [2026-02-09 23:xx] Startup dependency decoupling
- Problem
- Backend import path required `scanpy` at process startup, causing health check startup failure on lean environments.
- Root cause
- `scanpy` was imported at module import time in validator/runner service modules.
- Solution
- Converted `scanpy` imports to lazy runtime imports inside function scope.
- Code changes (files/functions)
- `backend/app/services/data_validator.py` (`validate_h5ad`)
- `backend/app/services/squidiff_runner.py` (`run_predict`)
- Verification results
- `python -c "from backend.app.main import app; print(app.title)"` prints app title successfully.
- `ruff`/`compileall` still pass.
- Impact assessment
- Service can boot earlier; validation/predict endpoints now return explicit dependency guidance if `scanpy` missing.

### [2026-02-10 00:xx] Frontend full workflow + Windows cmd_conda R support
- Problem
- Need to complete end-to-end UI flow and support user-specified Conda R environment on Windows CMD for Seurat conversion.
- Root cause
- Initial frontend was only a shell; backend R conversion only supported direct `Rscript`.
- Solution
- Implemented complete frontend workflow (upload/validate/train/poll/result).
- Added backend support for `cmd_conda` execution mode and per-request R runtime overrides.
- Added result asset URL mapping and job log API for result page display.
- Code changes (files/functions)
- `backend/app/core/config.py` (new `LABFLOW_R_*` settings)
- `backend/app/services/seurat_converter.py` (`_build_r_command`, `convert_to_h5ad`)
- `backend/app/api/datasets.py` (`ValidatePayload` extended with R runtime fields)
- `backend/app/api/jobs.py` (`GET /api/jobs/{job_id}/log`)
- `backend/app/api/results.py` (model detail API + asset URL + asset file serving)
- `frontend/src/services/api.ts` (full API client for flow)
- `frontend/src/App.tsx` (step-by-step workflow implementation)
- `frontend/src/styles/tokens.css` (form/status/result styles)
- `infra/.env.example` and `infra/docker-compose.yml` (runtime config exposure)
- Verification results
- Backend:
- `ruff check backend`: passed.
- `ruff format --check backend`: passed.
- `python -m compileall backend/app`: passed.
- Frontend:
- `npm install`: passed.
- `npm run lint`: passed.
- `npm run build`: passed.
- Impact assessment
- Lab users can now complete the requested workflow from a single frontend page.
- Windows R/Conda activation requirement is now configurable and executable via CMD.

### [2026-02-10 00:xx] V2 PRD drafted: Seurat interactive selection + 500x500 pipeline
- Problem
- Users still struggle with manual Seurat filtering and metadata preparation before training.
- Root cause
- MVP assumes preprocessed h5ad-style input and lacks in-UI cluster/group selection + bounded preprocessing.
- Solution
- Added a new PRD describing:
- 1) WebUI Seurat inspection and UMAP interaction,
- 2) user-selected metadata mapping (`group_column`, `cluster_column`),
- 3) cell stratified downsampling to max 500,
- 4) DEG-based gene selection to max 500,
- 5) final 500x500 training matrix contract.
- Code changes (files/functions)
- `docs/PRD_Seurat交互筛选与500x500训练管线.md` (new)
- Verification results
- Documentation update only (no runtime behavior changed in this step).
- Impact assessment
- V2 scope is now explicit and implementation-ready for phased development.

### [2026-02-10 01:xx] V2 Phase 1 initial implementation: Seurat inspect API + UI hook
- Problem
- Need to start implementing PRD V2 with API-first order, beginning from Seurat inspection capability.
- Root cause
- Existing MVP lacked a dedicated endpoint to expose metadata columns and UMAP preview from uploaded datasets.
- Solution
- Added backend Seurat inspection service and API endpoint: `POST /api/seurat/inspect`.
- Added frontend API client + new UI section to trigger inspection and render metadata/UMAP preview.
- Added API documentation and minimal backend endpoint tests.
- Code changes (files/functions)
- `backend/app/services/seurat_inspector.py` (`inspect_h5ad` and helpers)
- `backend/app/api/seurat.py` (`inspect_seurat` endpoint)
- `backend/app/main.py` (router registration)
- `frontend/src/services/api.ts` (`inspectSeurat`, inspect response types)
- `frontend/src/App.tsx` (new "Seurat 解析" step)
- `frontend/src/styles/tokens.css` (`chip-list`, `umap-preview`)
- `backend/tests/test_seurat_api.py` (404 + missing h5ad path checks)
- `docs/api/seurat.md` (new API doc)
- Verification results
- Frontend:
- `npm run lint`: passed.
- `npm run build`: passed.
- Backend smoke:
- `python` + `fastapi.testclient` script: `/api/health` 200 and `/api/seurat/inspect` missing dataset -> 404.
- Backend static check:
- `ruff check backend/app backend/tests`: passed (using custom `RUFF_CACHE_DIR` due default cache permission issue).
- Additional note:
- `uv run pytest ...` is still blocked by existing project dependency resolution constraints (`scanpy` + broad `requires-python`).
- Impact assessment
- PRD V2 Phase 1 now has a usable backend contract and frontend integration point.
- Full interactive筛选与500x500预处理（prepare-training）仍待后续 Phase 2/3 开发。

### [2026-02-10 02:xx] V2 Phase 2 implementation: prepare-training pipeline (500x500)
- Problem
- Need to implement PRD V2 Phase 2 backend pipeline with API contract: `/api/seurat/prepare-training` + job status query.
- Root cause
- Existing code can inspect Seurat metadata/UMAP but cannot produce bounded training matrix (`<=500 cells`, `<=500 genes`) with traceable reports.
- Solution
- Added dataset preprocessing service with:
- 1) cluster filtering by `selected_clusters`,
- 2) stratified sampling (`Group -> Cluster`) capped at 500 cells,
- 3) DEG-based gene selection (Wilcoxon) capped at 500 genes, with fallback to HVG/variance ranking.
- Added prepare-training APIs:
- `POST /api/seurat/prepare-training` for execution + dataset registration.
- `GET /api/seurat/prepare-training/{job_id}` for status/result.
- Added dedicated JSON state bucket for Seurat prepare jobs and API docs/tests.
- Code changes (files/functions)
- `backend/app/services/dataset_preprocessor.py`
- `stratified_sample_cells`, `select_top_genes`, `prepare_training_dataset`.
- `backend/app/api/seurat.py`
- `SeuratPrepareTrainingPayload`, `prepare_training`, `get_prepare_training_job`.
- `backend/app/storage/state_manager.py`
- new `seurat_prepare_jobs` store methods.
- `backend/tests/test_dataset_preprocessor.py`
- deterministic/bounded sampling tests.
- `backend/tests/test_seurat_api.py`
- prepare-training endpoint contract tests (error + success path with stubbed preprocessor).

### [2026-02-10 09:xx] WSL2 API simulation + dependency blockers
- Problem
- Need to run full API flow (upload/validate/inspect/prepare/train) using `data/coTC.rds` in WSL2.
- Root cause
- Backend could not boot without `python-multipart`, and sync FastAPI endpoints hung due to `anyio.to_thread.run_sync` in this environment.
- Solution
- Added multipart availability guard with `/api/datasets/register-local` for local-path testing.
- Converted API endpoints to `async def` to avoid anyio threadpool hangs under WSL2.
- Added explicit Rscript preflight check with clearer error message.
- Code changes (files/functions)
- `backend/app/api/datasets.py` (multipart guard + `register-local`; async endpoints)
- `backend/app/api/jobs.py` (async endpoints)
- `backend/app/api/results.py` (async endpoints)
- `backend/app/api/seurat.py` (async endpoints)
- `backend/app/main.py` (async health endpoint)
- `backend/app/services/seurat_converter.py` (Rscript availability check)
- `backend/app/services/visualize.py` (lazy imports for matplotlib/sklearn)
- `docs/api/datasets.md` (new API doc)
- Verification results
- ASGI smoke test via `httpx.AsyncClient`:
- `GET /api/health` -> 200 OK.
- `POST /api/datasets/register-local` -> 200 OK.
- `POST /api/datasets/{id}/validate` -> 400 with clear Rscript missing message.
- Checkfix:
- `ruff check backend/app backend/tests` -> passed.
- `ruff format --check backend/app backend/tests` -> passed (after formatting).
- Blockers
- `Rscript` not installed in WSL2, so `.rds` conversion fails.
- `scanpy` not installed, so `.h5ad` validation/inspect/prepare will fail after conversion.
- `python-multipart` cannot be installed due to offline pip, so multipart upload is disabled in this environment.
- Impact assessment
- API is usable in WSL2 for non-multipart endpoints; full Seurat pipeline requires R + SeuratDisk + scanpy installed.
- `docs/api/seurat.md`
- Phase 2 endpoints and payload/response docs.
- Verification results
- `ruff check backend/app backend/tests`: passed.
- `ruff format --check backend/app backend/tests`: passed.
- Backend smoke (with stubbed preprocessor): passed.
- `POST /api/seurat/prepare-training` returns `job_id` + `prepared_dataset_id`.
- `GET /api/seurat/prepare-training/{job_id}` returns `status=success`.
- Additional note
- `uv run pytest` remains blocked by existing dependency resolution issue (`scanpy` vs broad `requires-python` range), same as previous rounds.
- Impact assessment
- PRD V2 Phase 2 backend contract and core algorithm pipeline are now in place.
- Remaining PRD work is mainly Phase 3/4 (training flow默认接 prepared_dataset_id + frontend筛选页增强 + docs/UAT).

### [2026-02-10 02:xx] V2 Phase 3 implementation: train default prepared dataset + frontend summary
- Problem
- Need to make training jobs default to `prepared_dataset_id` and expose preprocessing source/summary in frontend.
- Root cause
- Train API previously always used incoming `dataset_id`, and frontend lacked prepare-summary context for training source traceability.
- Solution
- Backend train endpoint now resolves training dataset in priority:
- 1) if request `dataset_id` is already a prepared dataset, use itself;
- 2) else if `prepared_dataset_id` provided, validate it belongs to source dataset and use it;
- 3) else auto-pick latest prepared dataset derived from source dataset.
- Job metadata now includes `source_dataset_id`, `prepared_dataset_id`, `used_prepared_dataset`, plus param trace (`requested_dataset_id`, `train_dataset_id`).
- Frontend added prepare-training call/summary state and uses `prepared_dataset_id` by default when submitting train.
- Training and job status panels now show preprocessing summary and training source trace fields.
- Code changes (files/functions)
- `backend/app/api/jobs.py`
- `TrainJobPayload.prepared_dataset_id`, `_latest_prepared_dataset`, updated `submit_train_job`.
- `backend/tests/test_jobs_api.py`
- auto-select latest prepared + mismatched prepared id rejection tests.
- `frontend/src/services/api.ts`
- `prepareTraining` client, `PrepareTrainingResult`, extended `JobRecord`, train payload supports `preparedDatasetId`.
- `frontend/src/App.tsx`
- Phase 2 prepare form action, preprocessing summary display, default train-source wiring to prepared dataset, job trace fields.
- Verification results
- Backend checkfix:
- `ruff check backend/app backend/tests`: passed.
- `ruff format --check backend/app backend/tests`: passed.
- Frontend checkfix:
- `npm run lint`: passed.
- `npm run build`: passed.
- Backend smoke:
- train default prepared dataset selection script: passed (`phase3-train-default-smoke-ok`).
- Impact assessment
- Phase 3 core requirement is now met: train flow defaults to prepared dataset when available and source is visible in UI.
- Remaining items are mainly Phase 4 docs/UAT and richer交互筛选体验优化.

### [2026-02-10 03:xx] V2 Phase 4 implementation: docs completion + UAT delivery assets
- Problem
- Need to complete Phase 4 deliverables: V2 docs supplement, lab quickstart, and UAT script/checklist for at least two datasets.
- Root cause
- Existing docs covered base conversion and API but lacked a consolidated lab handoff package for V2 workflow and repeatable UAT execution.
- Solution
- Added V2 chapter to conversion guide (`docs/seurat转换指南.md`) with metadata规范、500x500约束、V2接口顺序与快速自检示例。
- Added lab handoff doc (`docs/实验室10分钟上手.md`) with practical timeline-oriented steps.
- Added executable UAT runner (`scripts/uat_phase4_seurat_v2.py`) supporting:
- repeated `--dataset-id` inputs (minimum two),
- inspect + prepare + optional train chain verification,
- bounded checks (`n_cells <= 500`, `n_genes <= 500`),
- JSON report output.
- Added checklist template (`docs/UAT_Seurat_V2_检查清单.md`) for manual acceptance tracking.
- Code changes (files/functions)
- `docs/seurat转换指南.md` (new V2 section)
- `docs/实验室10分钟上手.md` (new)
- `docs/UAT_Seurat_V2_检查清单.md` (new)
- `scripts/uat_phase4_seurat_v2.py` (`request_json`, `run_dataset_uat`, `poll_train_job_until_done`, CLI args)
- Verification results
- `python -m py_compile scripts/uat_phase4_seurat_v2.py`: passed.
- `ruff check scripts/uat_phase4_seurat_v2.py`: passed.
- `ruff format --check scripts/uat_phase4_seurat_v2.py`: passed.
- Impact assessment
- Phase 4 required delivery assets are now in repo and runnable.
- Lab members can run scripted UAT with two dataset IDs and archive JSON reports for handoff.

### [2026-02-10 03:xx] Documentation refresh: README + AGENTS + CLAUDE comprehensive rewrite
- Problem
- Need a complete, up-to-date project handbook covering front-end/back-end development, deployment, architecture, and feature status in one place.
- Root cause
- Existing top-level docs had partial overlap and outdated context (especially around V2 pipeline and local_Squidiff collaboration conventions).
- Solution
- Rewrote `README.md` as primary operator/developer entry:
- project capabilities, architecture, API map, local runbook, docker deployment, env vars, and doc index.
- Rewrote `AGENTS.md` as collaboration contract:
- layer boundaries, API-first workflow, module/API quick map, checkfix rules, doc sync rules, and skill trigger guidance.
- Rewrote `CLAUDE.md` as AI dev guide:
- current stack, backend/frontend architecture, V2 flow, quality gates, deployment and known constraints.
- Code changes (files/functions)
- `README.md` (full rewrite)
- `AGENTS.md` (full rewrite)
- `CLAUDE.md` (full rewrite)
- Verification results
- Documentation consistency scan passed:
- key API paths (`/api/seurat/prepare-training`, `/api/jobs/train`) and deployment commands are correctly reflected.
- Per user request, no runtime tests/checkfix were executed in this round (testing deferred to next day).
- Impact assessment
- New contributors and AI agents can now onboard with a single coherent set of docs for architecture, deployment and workflow expectations.

### [2026-02-11] Windows 全自动 500×500 测试 + cmd_conda / R 转换修复
- Problem
- 在 Windows 本机用 data/TC.rds（筋膜）、data/coTC.rds（结肠）模拟前端 API 跑 500×500 流程时，R 转换失败：conda.bat 在 cmd /c 下无法识别；转换进程返回 0 但未生成 h5ad。
- Root cause
- 1) cmd /c 单行命令中路径引号被解析成可执行名的一部分；2) SeuratDisk Convert() 生成的是 base.h5ad（替换 .h5seurat），R 脚本误用 paste0(..., ".h5ad") 得到 base.h5seurat.h5ad 导致 file.copy 失败；3) 传予 R 的 Windows 反斜杠路径在 R 中被转义，改用正斜杠可避免。
- Solution
- 1) cmd_conda 改为通过临时 .bat 文件执行（写入 call conda activate + Rscript ...），避免 cmd /c 引号问题；2) 向 R 传入路径时统一改为正斜杠；3) 临时 .bat 用毕删除；4) R 脚本中 converted 路径改为 sub("\\.h5seurat$", ".h5ad", tmp_h5seurat)，并对 file.copy 失败做 stop()；5) 新增 scripts/run_500x500_test_windows.py，模拟 register-local → validate（cmd_conda + r-4.3）→ inspect → prepare-training（自动推断 group/cluster 列），并校验 n_cells/n_genes ≤ 500。
- Code changes (files/functions)
- `backend/app/services/seurat_converter.py`（临时 .bat、正斜杠路径、错误信息含 stdout/stderr）
- `backend/scripts/seurat_to_h5ad.R`（converted 路径修正、file.copy 失败时 stop）
- `scripts/run_500x500_test_windows.py`（新建）
- `.debug/labflow-mvp-debug.md`（Windows 运行上下文、示例数据标注、测试脚本说明）
- Verification results
- 全自动测试：TC-筋膜、coTC-结肠 均完成 register → validate（R 转 h5ad）→ inspect → prepare-training，n_cells=500、n_genes=500，500×500 逻辑校验通过。
- Checkfix：`ruff check backend/app backend/tests` 通过；`ruff format backend/app` 已执行。
- Impact assessment
- Windows 下使用 conda r-4.3 的 R 转换与 500×500 流程可在本机一键脚本验证；后端需安装 scanpy 以支持 inspect。

### [2026-02-11] 全流程训练→预测→可视化报告脚本
- Problem
- 用户要求从 h5ad 转换开始到最终出可视化报告全部走通，metadata 由脚本/管线自行判断如何进入训练。
- Root cause
- 此前仅有 500×500 测试脚本，无训练、预测与报告拉取的一体化脚本。
- Solution
- 新增 `scripts/run_full_train_predict_viz_windows.py`：单条链路（默认 TC.rds）执行 register-local → validate（R 转 h5ad）→ inspect → prepare-training（500×500）→ POST /api/jobs/train（使用 prepared_dataset_id、gene_size/output_dim=500）→ 轮询训练完成 → POST /api/jobs/predict（同 prepared 数据 + 新 model_id）→ 轮询预测完成 → GET /api/results/job/{predict_job_id} → 下载 summary.assets 到 `scripts/output/full_flow_report/`（summary.json + pca_scatter.png、heatmap_top_var_genes.png）。metadata：prepared h5ad 的 .obs 已含 Group/Cluster，训练使用该 h5ad 表达矩阵（train_squidiff 不单独读 metadata）。
- Code changes (files/functions)
- `scripts/run_full_train_predict_viz_windows.py`（新建）
- `.gitignore`（增加 `scripts/output/`）
- Verification results
- 后端 `LABFLOW_DRY_RUN=true`、端口 8002 下执行脚本：转换→prepare→train→predict→报告下载 全部成功；报告目录含 summary.json 与 2 张 PNG。
- Checkfix：`ruff check` / `ruff format` 脚本通过。
- Impact assessment
- 一条命令可验证「转换→500×500→训练→预测→可视化报告」全流程；dry_run 下无需 GPU、数分钟内完成。

### [2026-02-11] 真实训练失败：rdkit 按需导入 + 训练错误信息增强
- Problem
- 关闭 LABFLOW_DRY_RUN 后跑全流程，训练子进程退出码 1，仅报 "Training command failed with exit code 1"，无法直接看到 train_squidiff.py 的报错。
- Root cause
- `Squidiff/scrna_datasets.py` 顶层 `from rdkit import Chem`，在 use_drug_structure=False 时也会触发导入，本机未安装 rdkit 导致 ModuleNotFoundError。Runner 未把子进程 stderr 写入异常信息。
- Solution
- 将 rdkit/Chem 导入移入 `Drug_dose_encoder` 内（仅 use_drug_structure=True 时调用），使无药物结构训练不依赖 rdkit。Runner 在训练失败时把 proc.stderr（或 stdout）末尾最多 1500 字写入 RuntimeError，便于 job error_msg 与日志排查。
- Code changes (files/functions)
- `Squidiff/scrna_datasets.py`（rdkit 按需导入）
- `backend/app/services/squidiff_runner.py`（run_train 失败时附带子进程输出）
- Verification results
- `ruff check` 通过。
- Impact assessment
- 无 rdkit 环境下 use_drug_structure=False 的训练可正常启动；后续若训练仍失败，错误信息会直接包含子进程输出，无需单独查 train.log。

### [2026-02-11] 训练轮询按本任务 PID 判断是否延长
- Problem
- 用户指出「nvidia-smi 进程列表中是否含 python 不够」，其它程序也可能用 GPU，应具体到**本脚本/本任务对应的训练进程**（进程号或路径）再判断是否延长。
- Root cause
- 脚本仅用「GPU 利用率 > 阈值」或「任意 Python 进程在 GPU」判断，易在多人/多任务共享 GPU 时误判。
- Solution
- 后端：`squidiff_runner.run_train` 改为用 `subprocess.Popen` 启动训练，得到子进程 PID 后通过新增参数 `on_start(pid)` 回调；`job_queue._execute_train` 在调用 run_train 时传入 `on_start=lambda pid: store.update_job(job_id, {"train_pid": pid})`，使 GET /api/jobs/{job_id} 返回 train_pid。脚本：新增 `get_gpu_pids()`（nvidia-smi --query-compute-apps=pid 及 -q 下 Process ID 解析），超时后先拉取 job；若存在 train_pid，则仅当 `train_pid in get_gpu_pids()` 时延长；否则退化为原逻辑（利用率或任意 Python 在 GPU）。延长/不延长时的提示改为「本任务训练进程 PID=xxx 仍在/已不在 GPU」。
- Code changes (files/functions)
- `backend/app/services/squidiff_runner.py`（run_train 增加 on_start，Popen + communicate 替代 run）
- `backend/app/services/job_queue.py`（_execute_train 传入 on_train_start 写 train_pid）
- `scripts/run_full_train_predict_viz_windows.py`（get_gpu_pids、轮询分支按 train_pid 判断延长）
- Verification results
- ruff check 通过；无新增 lint 报错。
- Impact assessment
- 仅在本任务训练进程（具体 PID）仍占用 GPU 时延长等待，其它程序占用 GPU 不会触发延长。

### [2026-02-11] Black 格式化与 Checkfix
- 用户已通过 `python -m black` 成功执行 black（PATH 已含 Python313\Scripts 或使用 python -m）。先对 `backend`、`scripts` 格式化 6 个文件；后对**全项目**执行 `python -m black .`，再格式化 Squidiff/ 与根目录共 14 个文件（dist_util.py、nn.py、respace.py、scrna_datasets.py、sample_squidiff.py、resample.py、Squidiff/train_squidiff.py、train_squidiff.py、MLPModel.py、fp16_util.py、logger.py、script_util.py、train_util.py、diffusion.py）。当前 `python -m black --check .` 与 `ruff check .` 均通过（45 files unchanged）。

### [2026-02-11] README 与部署文档精简（开发/部署只保留 uv、conda、python）
- 用户要求：开发部署只保留 uv、conda、直接 python，不整 uvicorn 等额外负担，和部署文档一起修干净。
- 修改：README 第 5 节改为「开发与部署（三种方式）」：5.1 环境准备仅列 uv / conda / 本机 Python 三选一；5.2 后端启动统一为 `python -m uvicorn backend.app.main:app --reload --host 0.0.0.0 --port 8000`；5.3 前端；5.4 Docker 一笔带过。部署文档：安装流程合并为 3.1 uv、3.2 conda、3.3 本机 Python、3.4 国内镜像；新增 5. 启动 LabFlow Web（一条后端命令 + 前端）；原 5–9 节顺延为 6–10。CLAUDE.md、AGENTS.md、backend/README.md 同步为「环境三选一 + python -m uvicorn」。

## Open Issues
- Real-world Seurat conversion relies on local R/SeuratDisk availability.
- Production auth is intentionally simplified for MVP.
- Full E2E Phase 2 run still depends on runtime `scanpy` + actual h5ad data availability.

## Technical Debt
- JSON file storage has limited concurrency compared with database-backed approach.
- Current UI is single-page workflow; multi-page routing and better UX states can be added later.
