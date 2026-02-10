# Datasets API

## POST /api/datasets/upload
Multipart upload for datasets (requires `python-multipart`).

Form fields
- `file` (required): dataset file (`.h5ad`, `.h5seurat`, `.rds`)
- `smiles_file` (optional): drug structure file
- `dataset_name` (optional): friendly name

Response
- `dataset`: dataset record

Notes
- If `python-multipart` is missing, this endpoint returns `503` with guidance.

## POST /api/datasets/register-local
Register a dataset from a local filesystem path (useful for local testing environments).

JSON body
- `local_path` (required): absolute or repo-relative path to dataset file
- `dataset_name` (optional)
- `smiles_path` (optional)
- `copy_to_upload_dir` (optional, default `true`): copy file into `LABFLOW_UPLOAD_DIR`

Response
- `dataset`: dataset record

## POST /api/datasets/{dataset_id}/validate
Validate (and optionally convert) the dataset into `.h5ad`.

JSON body
- `use_drug_structure` (optional, default `false`)
- `auto_convert` (optional, default `true`)
- `r_exec_mode` (optional)
- `r_conda_env` (optional)
- `r_conda_bat` (optional)
- `rscript_bin` (optional)

Response
- `dataset`
- `validation`

## GET /api/datasets
List datasets.

Response
- `items`: dataset records
