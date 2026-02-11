# README and Ruff Debug Record

## Metadata
- Module name: project-docs-maintenance
- Created: 2026-02-09 23:00:14 +08:00
- Last updated: 2026-02-09 23:10:13 +08:00
- Related files: `README.md`, `docs/部署文档.md`, `docs/seurat转换指南.md`, `docs/避坑指南.md`, `docs/环境安装指南.md`, `.debug/readme-ruff-debug.md`
- Dependency modules: training entry (`train_squidiff.py`), inference entry (`sample_squidiff.py`), CLI helper (`scripts/check_shape.py`)

## Runtime Context and Test Rules
- Runtime environment: local Windows terminal
- SSH mode (remote only): not used
- Remote project path (remote only): not used
- Verification and checkfix execution: run directly in local shell

## Context Relationship Network
- File structure:
  - `README.md` is the root user-facing documentation.
  - `train_squidiff.py` and `sample_squidiff.py` define the executable training and inference flow.
  - `scripts/check_shape.py` is the practical data-shape helper referenced by docs.
- Function call chain (doc-to-code mapping):
  - Training command -> `train_squidiff.py:parse_args()` -> `run_training()`
  - Inference sample -> `sample_squidiff.sampler(...)` -> `sampler.pred(...)`
- Variable and parameter dependency:
  - `gene_size` and `output_dim` must match data dimension.
  - `use_drug_structure=True` implies valid drug-related metadata and control dataset setup.
- Data flow:
  - `h5ad` input -> model encode/diffusion -> predicted expression output.

## Debug History
### [2026-02-11 16:45] README 用户视角 Quick Start 增补
- Problem statement:
  - User requested a user-perspective quick start section right after project introduction in `README.md`.
  - The section should clearly describe deployment and usage flow, and provide clickable doc links.
- Root cause:
  - Existing README had architecture and setup details, but no concise user-first "deploy then use" guide directly after introduction.
- Solution:
  - Added `## 1.5 Quick Start（用户视角）` immediately after section 1.
  - Structured flow into:
    - deployment path selection (local vs Docker)
    - startup sequence (backend/frontend and expected URLs)
    - end-user operation sequence (upload -> inspect -> prepare -> train -> result)
    - clickable docs map for next actions and troubleshooting
  - Added direct markdown links to:
    - `docs/部署文档.md`
    - `docs/seurat转换指南.md`
    - `docs/LabFlow前端用户操作说明.md`
    - `docs/实验室10分钟上手.md`
    - `docs/api/seurat.md`
    - `docs/api/jobs.md`
    - `docs/api/datasets.md`
    - `docs/UAT_Seurat_V2_检查清单.md`
    - `docs/模型能做什么与前端设计理念.md`
- Code changes:
  - `README.md`: inserted user-oriented quick-start section and doc links.
  - `.debug/readme-ruff-debug.md`: added this traceable record.
- Verification result:
  - Documentation-only update; no runtime code path changed.
  - IDE lint diagnostics were checked for edited files.
- Impact assessment:
  - Improves onboarding clarity for non-technical users.
  - Reduces time to first successful run in LabFlow UI.

### [2026-02-09 23:00] README cleanup and project-wide Ruff check
- Problem statement:
  - `README.md` had severe encoding corruption and inconsistent operational details.
  - User requested detailed README optimization plus whole-project Ruff validation.
- Root cause:
  - Existing documentation text included mojibake and mixed stale instructions.
  - `uv run ruff` path hit local cache permission issue, which blocked that route.
- Solution:
  - Rewrote `README.md` in clean ASCII with executable sections:
    - installation options
    - data requirements
    - quick-start training/inference
    - core training arguments
    - quality checks and project layout
  - Performed Ruff checkfix loop with fallback strategy:
    - Attempted `uv run ruff check .` -> failed with access error on uv cache.
    - Fallback `python -m ruff check .` -> passed.
- Code changes:
  - `README.md`: full rewrite for clarity and operability.
  - `.debug/readme-ruff-debug.md`: added this traceable record.
- Verification result:
  - `python -m ruff check .` -> `All checks passed!`
- Impact assessment:
  - No runtime Python code changed.
  - Documentation reliability improved for setup/train/inference workflows.
  - Lint status for the full project remains passing.

### [2026-02-09 23:10] Bilingual README rewrite and docs consolidation
- Problem statement:
  - User requested a higher-quality bilingual README with line-by-line EN/CN parallel text.
  - User requested docs optimization with merge of overlapping files and deletion of redundant docs.
- Root cause:
  - Root README was English-only and no longer met bilingual readability expectations.
  - `docs/部署文档.md` and `docs/环境安装指南.md` had substantial overlap.
  - Existing docs repeated similar setup and troubleshooting content across multiple files.
- Solution:
  - Rewrote `README.md` into EN/CN parallel style with aligned sections:
    - overview, capabilities, install, data requirements, quick start, core args, docs map
  - Consolidated docs:
    - merged install content into `docs/部署文档.md`
    - deleted redundant `docs/环境安装指南.md`
    - simplified `docs/seurat转换指南.md` to a practical conversion + validation flow
    - reduced `docs/避坑指南.md` to high-frequency troubleshooting and check templates
- Code changes:
  - `README.md`: bilingual parallel rewrite.
  - `docs/部署文档.md`: merged deployment + installation operations.
  - `docs/环境安装指南.md`: removed.
  - `docs/seurat转换指南.md`: simplified and de-duplicated.
  - `docs/避坑指南.md`: simplified and de-duplicated.
- Verification result:
  - `python -m ruff check .` -> `All checks passed!`
  - `python -m ruff format --check .` -> failed with pre-existing formatting drift in 15 Python files (outside this docs-only change scope)
- Impact assessment:
  - No runtime model/training logic changed.
  - Documentation structure is now shorter, less redundant, and easier to navigate.

## Pending Tracking Items
- Review whether additional docs under `docs/` should also be normalized to prevent future encoding confusion.

## Technical Debt Log
- `uv` execution path currently depends on local cache permissions; if required, clean cache policy or permission setup should be standardized.
- `ruff format --check` currently reports 15 files requiring formatting. This is pre-existing and should be handled in a dedicated formatting PR to avoid noisy mixed changes.

## Architecture Decision Notes
- This task is documentation-only for behavior clarity; no model architecture or training logic change introduced.
