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
- `docs/实验�?0分钟上手.md`
- `docs/UAT_Seurat_V2_检查清�?md`
- `scripts/uat_phase4_seurat_v2.py`
- `infra/docker-compose.yml`
- Dependency modules:
- `train_squidiff.py`
- `sample_squidiff.py`

## Runtime Context and Test Rules
- Runtime environment: **Windows 本机**（项目路�?`E:\Development\local_Squidiff`）；也可�?WSL2 下用 `/mnt/e/Development/local_Squidiff`�?
- SSH mode (if remote): Not used.
- Remote project path (if remote): N/A
- Validation/Checkfix execution mode: 在本�?PowerShell/CMD 直接执行；Windows �?R 须用 `cmd_conda` + `r-4.3`�?
- R execution constraint (confirmed by user): R conda env must be activated via `cmd` (not PowerShell). 推荐环境�?*r-4.3**（`F:\software\Miniconda3\envs\r-4.3`），包齐全、稳定�?
- R config strategy: support both `.env` defaults (`LABFLOW_R_*`) and per-request frontend overrides (`r_exec_mode`, `r_conda_env`, `r_conda_bat`, `rscript_bin`).
- 示例数据（data/）：**TC.rds** = 大鼠皮下筋膜针灸/痢疾 telocytes�?*coTC.rds** = 大鼠结肠针灸/痢疾 telocytes；细胞量较大，用�?500×500 流程测试�?
- Windows �?500×500 测试脚本：`scripts/run_500x500_test_windows.py`。先启动后端，再在另一终端执行 `python scripts/run_500x500_test_windows.py`；可选环境变�?`LABFLOW_BASE_URL`、`LABFLOW_R_CONDA_ENV`、`LABFLOW_R_CONDA_BAT`�?
- 全流程（转换→训练→预测→可视化报告）：`scripts/run_full_train_predict_viz_windows.py`。建议后端设�?`LABFLOW_DRY_RUN=true` 以快速跑通；报告输出�?`scripts/output/full_flow_report/`（summary.json + pca_scatter.png、heatmap_top_var_genes.png）�?
- **端口�?R 转换**：若 8000 被占用，可在其它端口启动后端并设 `LABFLOW_BASE_URL`。若 validate 报「conda.bat 不是内部或外部命令」，说明当前后端进程是旧代码，需**重启后端**以加�?R 转换的临�?.bat 修复（见 2026-02-11 Windows 全自�?500×500 测试条目）�?
- **前端校验�?Rscript was not found**：后端进�?PATH 中无 Rscript 时（�?R 仅在 Conda 环境 r-4.3 中），校验会 400。解决：在页面上将「R 执行方式」改�?**cmd_conda**，填�?conda.bat 完整路径（如 F:\\software\\Miniconda3\\condabin\\conda.bat）和 R 环境名（�?r-4.3）。已做：后端错误文案增加 cmd_conda 说明；前端校验区增加提示、并�?400 �?detail 解析后追加操作建议�?
- **前端用户说明�?*：已新增 `docs/LabFlow前端用户操作说明.md`，按页面 1) 上传 2) 校验 3) Seurat 解析 4) 500x500 预处�?5) 提交训练 6) 任务轮询 7) 结果�?逐项说明每个选项与参数如何填写（�?Windows Conda R、direct/cmd_conda、常见问题）。README �?5 节与�?9 节文档导航已引用该文档�?
- **前端即时反馈与运行日�?*：用户要求点按钮后要有反馈、长任务要有实时监测。已做：�?全局「当前任务」条带（indeterminate 进度�?+ 进行中文案），在 busyStep �?job �?queued/running 时显示；�?在「任务轮询状态」中增加「运行日志」小电视，轮�?GET /api/jobs/{id}/log �?2.5s，深色可滚动 pre 展示。详�?`.debug/ui-labflow-debug.md`�?
- **R/conda 参数菜单�?*：用户要求启动时探测 conda 环境，前端用菜单选择而非手输。已做：�?后端 GET /api/runtime/conda-envs（可�?query conda_bat）探�?Windows conda.bat 候选路径（PATH + 常见安装目录）并执行 conda env list 解析环境名；�?前端 2) 校验区：conda.bat 路径改为下拉（选项来自 API�?「其他（手动输入）」；Conda R 环境名改为下拉（选项来自 API）；页面加载时请�?conda-envs 并预填首候选与首个 R 风格环境（如 r-4.3）。校验失败时错误与推荐文案显示在按钮旁（step-error），不只在页顶�?
- **真实训练 vs dry_run**：后端设�?`LABFLOW_DRY_RUN=true` 时，训练不执行（只写占位 `model.pt`），预测用随机矩阵，图会正常生成。要得到真实训练出的模型，需**不设或关�?* `LABFLOW_DRY_RUN` 后重启后端再跑全流程；真实模型在 `backend/artifacts/jobs/<train_job_id>/checkpoints/` 下�?
- **训练失败 ModuleNotFoundError: rdkit**：当 `use_drug_structure=False` 时，Squidiff 不需 rdkit。已�?`Squidiff/scrna_datasets.py` 中将 rdkit 改为�?`Drug_dose_encoder` 内按需导入，避免无药物结构时因�?rdkit 导致训练启动失败。训练失败时后端会在 job �?error_msg �?train.log 中保留子进程 stderr 末尾，便于排查�?
- **训练 use_drug_structure 被误传为 True**：runner 原先�?`--use_drug_structure str(False)` �?`"False"`，argparse 解析�?True，导致脚本去读空�?control_data_path �?OSError。已改为仅当 `params["use_drug_structure"]` 为真时才追加 `--use_drug_structure True` �?`--control_data_path`，否则不传，使用 train_squidiff 默认 False。失败时若子进程�?stderr/stdout，则从已写入�?train.log 读末尾作为错误详情�?
- **训练轮询超时但显卡仍在跑**：脚本原为固�?3600s 超时，训练超�?1 小时即报错，�?GPU 仍在训练。已在全流程脚本中增加「超时后多侧面判断」：(1) nvidia-smi �?GPU 利用率，高于阈值则延长�?2) nvidia-smi 进程列表�?-query-compute-apps �?-q 解析）中若存在名称含 python 的进程，也视为训练可能仍在跑并延长。任一满足即延长等待（每次 30 分钟，总上�?4 小时）。环境变量：LABFLOW_TRAIN_GPU_BUSY_THRESHOLD、LABFLOW_TRAIN_EXTEND_SEC、LABFLOW_TRAIN_MAX_TOTAL_SEC�?
- **训练轮询按本任务 PID 判断（精确到进程�?*：不再仅看「是否有 python �?GPU」，改为优先�?*本任务训练进�?*是否仍在 GPU。后端在启动训练子进程时�?Popen 获得 PID，通过 on_start(pid) 回调写入 job �?train_pid；GET /api/jobs/{job_id} 返回该字段。脚本增�?get_gpu_pids()（nvidia-smi --query-compute-apps=pid �?-q 解析 Process ID），超时后若 job �?train_pid，则**仅当 train_pid in get_gpu_pids()** 时才延长；否则退化为「利用率或任�?Python �?GPU」。这样其它程序占�?GPU 不会误触发延长�?

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
- Full interactive筛选与500x500预处理（prepare-training）仍待后�?Phase 2/3 开发�?

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
- Remaining PRD work is mainly Phase 3/4 (training flow默认�?prepared_dataset_id + frontend筛选页增强 + docs/UAT).

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
- Remaining items are mainly Phase 4 docs/UAT and richer交互筛选体验优�?

### [2026-02-10 03:xx] V2 Phase 4 implementation: docs completion + UAT delivery assets
- Problem
- Need to complete Phase 4 deliverables: V2 docs supplement, lab quickstart, and UAT script/checklist for at least two datasets.
- Root cause
- Existing docs covered base conversion and API but lacked a consolidated lab handoff package for V2 workflow and repeatable UAT execution.
- Solution
- Added V2 chapter to conversion guide (`docs/seurat转换指南.md`) with metadata规范�?00x500约束、V2接口顺序与快速自检示例�?
- Added lab handoff doc (`docs/实验�?0分钟上手.md`) with practical timeline-oriented steps.
- Added executable UAT runner (`scripts/uat_phase4_seurat_v2.py`) supporting:
- repeated `--dataset-id` inputs (minimum two),
- inspect + prepare + optional train chain verification,
- bounded checks (`n_cells <= 500`, `n_genes <= 500`),
- JSON report output.
- Added checklist template (`docs/UAT_Seurat_V2_检查清�?md`) for manual acceptance tracking.
- Code changes (files/functions)
- `docs/seurat转换指南.md` (new V2 section)
- `docs/实验�?0分钟上手.md` (new)
- `docs/UAT_Seurat_V2_检查清�?md` (new)
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

### [2026-02-11] Windows 全自�?500×500 测试 + cmd_conda / R 转换修复
- Problem
- �?Windows 本机�?data/TC.rds（筋膜）、data/coTC.rds（结肠）模拟前端 API �?500×500 流程时，R 转换失败：conda.bat �?cmd /c 下无法识别；转换进程返回 0 但未生成 h5ad�?
- Root cause
- 1) cmd /c 单行命令中路径引号被解析成可执行名的一部分�?) SeuratDisk Convert() 生成的是 base.h5ad（替�?.h5seurat），R 脚本误用 paste0(..., ".h5ad") 得到 base.h5seurat.h5ad 导致 file.copy 失败�?) 传予 R �?Windows 反斜杠路径在 R 中被转义，改用正斜杠可避免�?
- Solution
- 1) cmd_conda 改为通过临时 .bat 文件执行（写�?call conda activate + Rscript ...），避免 cmd /c 引号问题�?) �?R 传入路径时统一改为正斜杠；3) 临时 .bat 用毕删除�?) R 脚本�?converted 路径改为 sub("\\.h5seurat$", ".h5ad", tmp_h5seurat)，并�?file.copy 失败�?stop()�?) 新增 scripts/run_500x500_test_windows.py，模�?register-local �?validate（cmd_conda + r-4.3）→ inspect �?prepare-training（自动推�?group/cluster 列），并校验 n_cells/n_genes �?500�?
- Code changes (files/functions)
- `backend/app/services/seurat_converter.py`（临�?.bat、正斜杠路径、错误信息含 stdout/stderr�?
- `backend/scripts/seurat_to_h5ad.R`（converted 路径修正、file.copy 失败�?stop�?
- `scripts/run_500x500_test_windows.py`（新建）
- `.debug/labflow-mvp-debug.md`（Windows 运行上下文、示例数据标注、测试脚本说明）
- Verification results
- 全自动测试：TC-筋膜、coTC-结肠 均完�?register �?validate（R �?h5ad）→ inspect �?prepare-training，n_cells=500、n_genes=500�?00×500 逻辑校验通过�?
- Checkfix：`ruff check backend/app backend/tests` 通过；`ruff format backend/app` 已执行�?
- Impact assessment
- Windows 下使�?conda r-4.3 �?R 转换�?500×500 流程可在本机一键脚本验证；后端需安装 scanpy 以支�?inspect�?

### [2026-02-11] 全流程训练→预测→可视化报告脚本
- Problem
- 用户要求�?h5ad 转换开始到最终出可视化报告全部走通，metadata 由脚�?管线自行判断如何进入训练�?
- Root cause
- 此前仅有 500×500 测试脚本，无训练、预测与报告拉取的一体化脚本�?
- Solution
- 新增 `scripts/run_full_train_predict_viz_windows.py`：单条链路（默认 TC.rds）执�?register-local �?validate（R �?h5ad）→ inspect �?prepare-training�?00×500）→ POST /api/jobs/train（使�?prepared_dataset_id、gene_size/output_dim=500）→ 轮询训练完成 �?POST /api/jobs/predict（同 prepared 数据 + �?model_id）→ 轮询预测完成 �?GET /api/results/job/{predict_job_id} �?下载 summary.assets �?`scripts/output/full_flow_report/`（summary.json + pca_scatter.png、heatmap_top_var_genes.png）。metadata：prepared h5ad �?.obs 已含 Group/Cluster，训练使用该 h5ad 表达矩阵（train_squidiff 不单独读 metadata）�?
- Code changes (files/functions)
- `scripts/run_full_train_predict_viz_windows.py`（新建）
- `.gitignore`（增�?`scripts/output/`�?
- Verification results
- 后端 `LABFLOW_DRY_RUN=true`、端�?8002 下执行脚本：转换→prepare→train→predict→报告下�?全部成功；报告目录含 summary.json �?2 �?PNG�?
- Checkfix：`ruff check` / `ruff format` 脚本通过�?
- Impact assessment
- 一条命令可验证「转换→500×500→训练→预测→可视化报告」全流程；dry_run 下无需 GPU、数分钟内完成�?

### [2026-02-11] 真实训练失败：rdkit 按需导入 + 训练错误信息增强
- Problem
- 关闭 LABFLOW_DRY_RUN 后跑全流程，训练子进程退出码 1，仅�?"Training command failed with exit code 1"，无法直接看�?train_squidiff.py 的报错�?
- Root cause
- `Squidiff/scrna_datasets.py` 顶层 `from rdkit import Chem`，在 use_drug_structure=False 时也会触发导入，本机未安�?rdkit 导致 ModuleNotFoundError。Runner 未把子进�?stderr 写入异常信息�?
- Solution
- �?rdkit/Chem 导入移入 `Drug_dose_encoder` 内（�?use_drug_structure=True 时调用），使无药物结构训练不依赖 rdkit。Runner 在训练失败时�?proc.stderr（或 stdout）末尾最�?1500 字写�?RuntimeError，便�?job error_msg 与日志排查�?
- Code changes (files/functions)
- `Squidiff/scrna_datasets.py`（rdkit 按需导入�?
- `backend/app/services/squidiff_runner.py`（run_train 失败时附带子进程输出�?
- Verification results
- `ruff check` 通过�?
- Impact assessment
- �?rdkit 环境�?use_drug_structure=False 的训练可正常启动；后续若训练仍失败，错误信息会直接包含子进程输出，无需单独�?train.log�?

### [2026-02-11] 训练轮询按本任务 PID 判断是否延长
- Problem
- 用户指出「nvidia-smi 进程列表中是否含 python 不够」，其它程序也可能用 GPU，应具体�?*本脚�?本任务对应的训练进程**（进程号或路径）再判断是否延长�?
- Root cause
- 脚本仅用「GPU 利用�?> 阈值」或「任�?Python 进程�?GPU」判断，易在多人/多任务共�?GPU 时误判�?
- Solution
- 后端：`squidiff_runner.run_train` 改为�?`subprocess.Popen` 启动训练，得到子进程 PID 后通过新增参数 `on_start(pid)` 回调；`job_queue._execute_train` 在调�?run_train 时传�?`on_start=lambda pid: store.update_job(job_id, {"train_pid": pid})`，使 GET /api/jobs/{job_id} 返回 train_pid。脚本：新增 `get_gpu_pids()`（nvidia-smi --query-compute-apps=pid �?-q �?Process ID 解析），超时后先拉取 job；若存在 train_pid，则仅当 `train_pid in get_gpu_pids()` 时延长；否则退化为原逻辑（利用率或任�?Python �?GPU）。延�?不延长时的提示改为「本任务训练进程 PID=xxx 仍在/已不�?GPU」�?
- Code changes (files/functions)
- `backend/app/services/squidiff_runner.py`（run_train 增加 on_start，Popen + communicate 替代 run�?
- `backend/app/services/job_queue.py`（_execute_train 传入 on_train_start �?train_pid�?
- `scripts/run_full_train_predict_viz_windows.py`（get_gpu_pids、轮询分支按 train_pid 判断延长�?
- Verification results
- ruff check 通过；无新增 lint 报错�?
- Impact assessment
- 仅在本任务训练进程（具体 PID）仍占用 GPU 时延长等待，其它程序占用 GPU 不会触发延长�?

### [2026-02-11] Black 格式化与 Checkfix
- 用户已通过 `python -m black` 成功执行 black（PATH 已含 Python313\Scripts 或使�?python -m）。先�?`backend`、`scripts` 格式�?6 个文件；后对**全项�?*执行 `python -m black .`，再格式�?Squidiff/ 与根目录�?14 个文件（dist_util.py、nn.py、respace.py、scrna_datasets.py、sample_squidiff.py、resample.py、Squidiff/train_squidiff.py、train_squidiff.py、MLPModel.py、fp16_util.py、logger.py、script_util.py、train_util.py、diffusion.py）。当�?`python -m black --check .` �?`ruff check .` 均通过�?5 files unchanged）�?

### [2026-02-11] README 与部署文档精简（开�?部署只保�?uv、conda、python�?
- 用户要求：开发部署只保留 uv、conda、直�?python，不�?uvicorn 等额外负担，和部署文档一起修干净�?
- 修改：README �?5 节改为「开发与部署（三种方式）」：5.1 环境准备仅列 uv / conda / 本机 Python 三选一�?.2 后端启动统一�?`python -m uvicorn backend.app.main:app --reload --host 0.0.0.0 --port 8000`�?.3 前端�?.4 Docker 一笔带过。部署文档：安装流程合并�?3.1 uv�?.2 conda�?.3 本机 Python�?.4 国内镜像；新�?5. 启动 LabFlow Web（一条后端命�?+ 前端）；�?5�? 节顺延为 6�?0。CLAUDE.md、AGENTS.md、backend/README.md 同步为「环境三选一 + python -m uvicorn」�?

## Open Issues
- Real-world Seurat conversion relies on local R/SeuratDisk availability.
- Production auth is intentionally simplified for MVP.
- Full E2E Phase 2 run still depends on runtime `scanpy` + actual h5ad data availability.

## Technical Debt
- JSON file storage has limited concurrency compared with database-backed approach.
- Current UI is single-page workflow; multi-page routing and better UX states can be added later.

### [2026-02-11 19:16] Docs fix: install CUDA torch before requirements
- Problem
- Backend env setup could install CPU-only torch first when running requirements before explicit CUDA torch install.
- Root cause
- Install commands in README/deployment docs were either missing explicit torch step or had the order reversed.
- Solution
- Updated install flow for uv/conda/venv to install CUDA torch first, then install requirements with --extra-index-url pointing to PyTorch CUDA wheels to avoid fallback/override to CPU builds.
- Code changes (files/functions)
- README.md (section 5.1 env setup)
- docs/�����ĵ�.md (sections 3.1/3.2/3.3)
- Verification results
- Manual doc diff check passed; command order is now CUDA torch first in both docs.
- Impact assessment
- Reduces risk of accidental CPU-only PyTorch installation in backend setup; keeps deploy guidance consistent across docs.


### [2026-02-11 19:32] Per-user scheduler mode (serial/parallel) with default parallel cap=3
- Problem
- Existing queue model was effectively single-thread serial; users needed a simple per-user switch between serial and parallel behavior without admin console complexity.
- Root cause
- JobQueue had a single worker and no user-level concurrency policy; no user preference API existed.
- Solution
- Added per-user scheduler preference API (`/api/user-prefs/scheduler`) backed by SQLite (`user_prefs` table), stored `owner_user_id` on submitted jobs, and refactored JobQueue into multi-worker scheduling with per-user limits: serial=1, parallel=3.
- Code changes (files/functions)
- `backend/app/services/auth_service.py` (user_prefs schema + get/set scheduler mode)
- `backend/app/api/user_prefs.py` (new API)
- `backend/app/api/jobs.py` (persist owner_user_id)
- `backend/app/services/job_queue.py` (multi-worker dispatcher + per-user slot control)
- `backend/app/main.py` + `backend/app/runtime.py` (router/runtime wiring)
- `frontend/src/services/api.ts` + `frontend/src/App.tsx` (scheduler mode UI + API integration)
- `docs/api/jobs.md` + `docs/api/user_prefs.md` + `docs/LabFlowǰ���û�����˵��.md` + `docs/�����ĵ�.md` + `README.md`
- `backend/tests/test_user_prefs_api.py` (new tests)
- Verification results
- `python -m compileall backend/app` passed.
- `python -m pytest ...` blocked in current environment (`ModuleNotFoundError: fastapi`).
- `ruff` blocked (ruff not installed in current environment).
- `npm run lint` executed; `npm run build` blocked by host permission issue (`spawn EPERM`).
- Impact assessment
- Each logged-in user can now choose task scheduling mode in Task Center: serial(1) or parallel(3). Queue dispatch is user-isolated by owner_user_id; changing mode affects scheduling immediately for queued tasks.

### [2026-02-11 19:53] Auth guide unreachable + login click appears unresponsive (network/error clarity hardening)
- Problem
- User reported /api/auth/user-guide�򲻿� and login actions seemed to have no response in frontend.
- Root cause
- Two risk points were identified: (1) guide file discovery could hit non-target markdown and fail hard in render path; (2) frontend fetch lacked timeout and explicit backend-unreachable messaging, so network stalls looked like no-op.
- Solution
- Backend: guide path now prefers docs/LabFlowǰ���û�����˵��.md; render failure falls back to raw markdown response.
- Frontend: added request timeout (15s) and normalized network error message (Cannot reach backend API ...) for auth/upload/general API calls.
- Code changes (files/functions)
- ackend/app/api/auth.py (_find_user_guide_path, user_guide)
- rontend/src/services/api.ts (etchWithTimeout, 	oNetworkError, requestJson/uploadDataset callsites)
- Verification results
- 
pm run lint executed.
- Direct network probe to http://192.168.1.104:8000/* from current environment timed out (indicates backend/network reachability issue in this execution context).
- Impact assessment
- User guide endpoint is more robust; frontend now fails fast with explicit backend/network guidance instead of silent waiting behavior.


### [2026-02-11 20:06] README �����޸�������ͳһΪ UTF-8��
- Problem
- README �ڲ��ֻ������� GitHub/�༭������ʾ���롣
- Root cause
- README �ļ����벻�� UTF-8�����¿�ƽ̨��ȡʱ���벻һ�¡�
- Solution
- �� README.md �ӱ������ı�������ת��Ϊ UTF-8��LF����
- Code changes (files/functions)
- README.md������ת�������������ݸĶ���
- Verification results
- Python �ϸ�У�飺README.md �ɱ� utf-8 �ɹ����롣
- �ؼ���������У��ͨ������ CJK �ַ�����
- Impact assessment
- README �� GitHub �� UTF-8 �������¿��ȶ���ʾ���ģ����ٳ������롣


### [2026-02-11 20:09] README ��дΪ��Ӣ˫�����һ�� + Ӣ��һ�У�
- Problem
- �û�Ҫ�� README ������Ӣ˫���Ű棨����һ�С�Ӣ��һ�У������޸���ʷ����Ӱ�졣
- Root cause
- README ֮ǰ���ڱ�����������ݽṹ��ͳһ�����¿ɶ��Բ
- Solution
- ֱ����д README.md Ϊͳһ˫��ṹ������ logo��ԭ��Ŀ��л����װ�������ĵ�������
- Code changes (files/functions)
- README.md��ȫ����д��UTF-8��
- Verification results
- Python ��ȡУ��ͨ����README.md �ɰ� UTF-8 ���룬�������ݴ��ڡ�
- Impact assessment
- README ������/Ӣ�Ķ��߸��Ѻã��ұ���绷���������⡣


### [2026-02-11 20:35] Launcher switched to uv-run backend + Windows packaging command fix
- Problem
- Running `python labflow_launcher.py` started backend with `python -m uvicorn ...` from Anaconda base, which failed with `No module named uvicorn`.
- Running `pyinstaller ...` also failed because `pyinstaller.exe` was installed in user Scripts but not on PATH.
- Root cause
- Launcher backend command was hardcoded to interpreter-local uvicorn module (`python -m uvicorn`) instead of uv-managed execution.
- Packaging docs used `pyinstaller` executable directly, which is PATH-sensitive on Windows.
- Solution
- Changed launcher backend startup command to `uv run python -m uvicorn backend.app.main:app ...`.
- Added `LABFLOW_UV` support to allow full-path uv binary when uv is not on PATH.
- Updated docs to use `python -m PyInstaller ...` for stable packaging on Windows.
- Updated README and deployment docs backend start command to uv-run style for consistency.
- Code changes (files/functions)
- `labflow_launcher.py` (`LauncherConfig`, `detect_uv_command`, `build_backend_cmd`, startup config in `main`)
- `README.md` (backend startup command)
- `docs/�����ĵ�.md` (backend startup command + auth example)
- `docs/Windowsһ��������.md` (backend command, packaging command, troubleshooting, LABFLOW_UV env)
- Verification results
- `python -m py_compile labflow_launcher.py`: passed.
- `python labflow_launcher.py --dry-run`: now fails fast with clear uv guidance in current shell (`Missing command uv ... set LABFLOW_UV`).
- `python -m ruff ...`: blocked in current environment (`No module named ruff`).
- Impact assessment
- Launcher no longer depends on base-conda uvicorn availability and follows uv-first workflow.
- Windows packaging instructions no longer depend on PATH containing `pyinstaller.exe`.
- Remaining prerequisite: uv must be available via PATH or `LABFLOW_UV`.
### [2026-02-11 20:49] uv resolver split failure on launcher startup (scanpy vs requires-python>=3.8)
- Problem
- User reran launcher and backend startup failed before health check. uv reported unsatisfiable resolution across Python split markers, driven by `scanpy>=1.10.0` and project `requires-python >=3.8`.
- Root cause
- `uv run` in project mode resolves project metadata across declared Python range. Current range included Python 3.8, but scanpy in requested range does not support 3.8, causing resolver failure before command execution.
- Solution
- Launcher backend command switched to `uv run --active --no-project ...` so it uses the already activated environment and skips project metadata resolution for startup.
- Tightened project metadata to `requires-python = ">=3.9"` to align with dependency compatibility and avoid future resolver conflicts.
- Synced startup docs to the same `uv run --active --no-project` command.
- Code changes (files/functions)
- `labflow_launcher.py` (`build_backend_cmd`)
- `pyproject.toml` (`requires-python`)
- `README.md` (backend startup command)
- `docs/�����ĵ�.md` (backend startup command and auth example)
- `docs/Windowsһ��������.md` (launcher backend command)
- Verification results
- `python -m py_compile labflow_launcher.py`: passed.
- String checks confirm backend command now includes `--active --no-project` in code and docs.
- ruff check unavailable in this shell (`No module named ruff`).
- Impact assessment
- Launcher startup no longer blocks on uv project dependency resolution splits.
- Project metadata is now consistent with scanpy's minimum supported Python range.
### [2026-02-11 21:02] Frontend spawn WinError 2 on launcher (npm command resolution hardening)
- Problem
- Backend started and passed health check, but launcher failed immediately when spawning frontend with `[WinError 2] ϵͳ�Ҳ���ָ�����ļ�`.
- Root cause
- Frontend startup used bare `npm` command token. In some Windows shells/environments this can pass pre-check but still fail on subprocess spawn due command resolution differences.
- Solution
- Added explicit npm command resolution via `detect_npm_command()`.
- Launcher now stores and uses an absolute npm executable path (`npm.cmd`) for both frontend build and frontend run.
- Added `LABFLOW_NPM` env override for environments where npm resolution is unstable.
- Updated Windows launcher doc with `LABFLOW_NPM` and WinError 2 troubleshooting.
- Code changes (files/functions)
- `labflow_launcher.py` (`LauncherConfig.npm_cmd`, `detect_npm_command`, `prepare_frontend_if_needed`, `build_frontend_cmd`, `main` config assembly)
- `docs/Windowsһ��������.md` (env var and troubleshooting entries)
- Verification results
- `python -m py_compile labflow_launcher.py`: passed.
- Source checks confirm frontend commands now use `config.npm_cmd`.
- Impact assessment
- Reduces Windows command resolution failures when starting frontend from Python subprocess.
### [2026-02-11 21:14] Frontend startup failed: missing npm preview script
- Problem
- Launcher reached frontend spawn stage, but npm returned `Missing script: "preview"`, causing launcher shutdown.
- Root cause
- `frontend/package.json` had `dev` and `build`, but no `preview` script, while launcher default mode is `preview`.
- Solution
- Added `"preview": "vite preview"` script to frontend package scripts.
- Code changes (files/functions)
- `frontend/package.json` (`scripts.preview`)
- Verification results
- In this execution sandbox, direct frontend build/preview verification is constrained (EPERM on esbuild spawn), but script registration is confirmed and launcher will no longer fail with `Missing script: preview`.
- Impact assessment
- Launcher default frontend mode (`preview`) now matches frontend scripts and can proceed in normal Windows runtime environments.

### [2026-02-12 16:xx] Validate step false backend unreachable message
- Problem
- Frontend validate step showed `Cannot reach backend API ...` while browser network showed `POST /api/datasets/{id}/validate` returned HTTP 400.
- Root cause
- Request timeout and network failures shared a single generic message in `frontend/src/services/api.ts`.
- Validate request can run longer (R/Conda conversion), so timeout was misreported as backend unreachable.
- Solution
- Added per-request timeout support to request client.
- Split timeout error message from network-unreachable message.
- Increased dataset validate timeout to 180s.
- Code changes (files/functions)
- `frontend/src/services/api.ts`: `fetchWithTimeout`, `requestJson`, `toNetworkError`, `validateDataset`.
- Verification results
- Frontend checkfix in progress in this round (`npm run lint`, `npm run build`).
- Impact assessment
- Validate step now gives accurate timeout feedback and is less likely to fail early during slow R conversion.

### [2026-02-12 17:xx] Large dataset inspect stalls: timeout + memory pressure mitigation
- Problem
- User reported large datasets often stop at validate/Seurat inspect and cannot continue.
- Root cause
- Frontend request timeout was shorter than heavy processing in some paths.
- Seurat inspect loaded h5ad in full memory mode, which is fragile for large files.
- Solution
- Frontend: unified API timeout to 10 minutes (600s) for all requests.
- Backend: seurat inspector now reads h5ad in backed mode (`backed="r"`) and closes file handle after inspect.
- Code changes (files/functions)
- `frontend/src/services/api.ts`: `REQUEST_TIMEOUT_MS=600000`, `VALIDATE_TIMEOUT_MS=600000`.
- `backend/app/services/seurat_inspector.py`: `_load_adata`, `_close_adata`, `inspect_h5ad` resource lifecycle.
- `docs/LabFlowǰ���û�����˵��.md`: added troubleshooting section for large datasets.
- Verification results
- `npm run lint`: passed in this environment.
- `npm run build`: blocked by environment permission (`spawn EPERM` from esbuild).
- `ruff check backend/app backend/tests`: command unavailable in current environment.
- `ruff format --check backend/app backend/tests`: command unavailable in current environment.
- Impact assessment
- Long-running validate/inspect requests are less likely to fail early due to frontend timeout.
- Large h5ad inspect memory risk is reduced by backed read mode.
