export type HealthResponse = {
  status: string;
};

export type ValidationReport = {
  valid: boolean;
  errors: string[];
  warnings: string[];
  summary?: {
    n_cells: number;
    n_genes: number;
    obs_columns: string[];
  };
};

export type UmapPoint = {
  cell_id: string;
  x: number;
  y: number;
};

export type SeuratInspectReport = {
  n_cells: number;
  n_genes: number;
  metadata_columns: string[];
  /** 每列的唯一取值及计数（类似 R table()），便于点击标签查看分类（如 celltype 下有哪些类型及细胞数） */
  metadata_column_values?: Record<
    string,
    Array<{ value: string; count: number }>
  >;
  has_umap: boolean;
  umap?: {
    key: string;
    n_points: number;
    preview_count: number;
    truncated: boolean;
    preview: UmapPoint[];
  } | null;
  warnings: string[];
};

export type PrepareTrainingResult = {
  job_id: string;
  prepared_dataset_id: string;
  n_cells: number;
  n_genes: number;
  sampling_report: {
    mode: string;
    seed: number;
    input_cells: number;
    output_cells: number;
    [key: string]: unknown;
  };
  gene_report: {
    method: string;
    input_genes: number;
    output_genes: number;
    fallback_used?: boolean;
    fallback_reason?: string | null;
    selected_genes?: string[];
    [key: string]: unknown;
  };
};

export type DatasetRecord = {
  id: string;
  name: string;
  status: string;
  input_type: string;
  path_raw: string;
  path_h5ad?: string | null;
  smiles_path?: string | null;
  validation?: ValidationReport | null;
  created_at: string;
  updated_at: string;
};

export type JobRecord = {
  id: string;
  type: "train" | "predict";
  status: "queued" | "running" | "success" | "failed" | "canceled";
  dataset_id: string;
  source_dataset_id?: string;
  prepared_dataset_id?: string | null;
  used_prepared_dataset?: boolean;
  params: Record<string, unknown>;
  model_id?: string;
  result_id?: string;
  error_msg?: string | null;
  log_path?: string;
  started_at?: string | null;
  ended_at?: string | null;
  created_at: string;
  updated_at: string;
};

export type ModelRecord = {
  id: string;
  job_id: string;
  dataset_id: string;
  path_ckpt: string;
  log_path: string;
  gene_size: number;
  output_dim: number;
  use_drug_structure: boolean;
  metrics: Record<string, unknown>;
  created_at: string;
  updated_at: string;
};

export type ResultAsset = {
  type: string;
  path: string;
  name: string;
  url: string;
};

export type ResultRecord = {
  id: string;
  job_id: string;
  dataset_id: string;
  prediction_path: string;
  log_path: string;
  artifact_dir: string;
  summary: {
    n_cells: number;
    n_genes: number;
    mean_expression: number;
    std_expression: number;
    assets: ResultAsset[];
  };
  created_at: string;
  updated_at: string;
};

const API_BASE = import.meta.env.VITE_API_BASE ?? "http://localhost:8000";

type JsonMethod = "GET" | "POST";

async function requestJson<T>(
  path: string,
  method: JsonMethod = "GET",
  body?: unknown
): Promise<T> {
  const response = await fetch(`${API_BASE}${path}`, {
    method,
    headers: body ? { "Content-Type": "application/json" } : undefined,
    body: body ? JSON.stringify(body) : undefined
  });

  if (!response.ok) {
    const text = await response.text();
    throw new Error(text || `${method} ${path} failed: ${response.status}`);
  }
  return (await response.json()) as T;
}

export function buildApiUrl(path: string): string {
  return `${API_BASE}${path}`;
}

export async function getHealth(): Promise<HealthResponse> {
  return requestJson<HealthResponse>("/api/health");
}

export async function uploadDataset(input: {
  datasetFile: File;
  datasetName?: string;
  smilesFile?: File | null;
}): Promise<DatasetRecord> {
  const formData = new FormData();
  formData.append("file", input.datasetFile);
  if (input.datasetName) {
    formData.append("dataset_name", input.datasetName);
  }
  if (input.smilesFile) {
    formData.append("smiles_file", input.smilesFile);
  }
  const response = await fetch(`${API_BASE}/api/datasets/upload`, {
    method: "POST",
    body: formData
  });
  if (!response.ok) {
    throw new Error(await response.text());
  }
  const payload = (await response.json()) as { dataset: DatasetRecord };
  return payload.dataset;
}

export async function validateDataset(input: {
  datasetId: string;
  useDrugStructure: boolean;
  autoConvert?: boolean;
  rExecMode?: "direct" | "cmd_conda";
  rCondaEnv?: string;
  rCondaBat?: string;
  rscriptBin?: string;
}): Promise<{ dataset: DatasetRecord; validation: ValidationReport }> {
  return requestJson(`/api/datasets/${input.datasetId}/validate`, "POST", {
    use_drug_structure: input.useDrugStructure,
    auto_convert: input.autoConvert ?? true,
    r_exec_mode: input.rExecMode,
    r_conda_env: input.rCondaEnv,
    r_conda_bat: input.rCondaBat,
    rscript_bin: input.rscriptBin
  });
}

export async function inspectSeurat(input: {
  datasetId: string;
  umapPreviewLimit?: number;
}): Promise<SeuratInspectReport> {
  const payload = await requestJson<{
    dataset_id: string;
    inspect: SeuratInspectReport;
  }>("/api/seurat/inspect", "POST", {
    dataset_id: input.datasetId,
    umap_preview_limit: input.umapPreviewLimit ?? 500
  });
  return payload.inspect;
}

export async function prepareTraining(input: {
  datasetId: string;
  groupColumn: string;
  clusterColumn: string;
  selectedClusters: string[];
  seed: number;
}): Promise<PrepareTrainingResult> {
  return requestJson<PrepareTrainingResult>("/api/seurat/prepare-training", "POST", {
    dataset_id: input.datasetId,
    group_column: input.groupColumn,
    cluster_column: input.clusterColumn,
    selected_clusters: input.selectedClusters,
    seed: input.seed
  });
}

export async function createTrainJob(input: {
  datasetId: string;
  preparedDatasetId?: string;
  geneSize: number;
  outputDim: number;
  useDrugStructure: boolean;
  batchSize: number;
  lr: number;
}): Promise<JobRecord> {
  const payload = await requestJson<{ job: JobRecord }>("/api/jobs/train", "POST", {
    dataset_id: input.datasetId,
    prepared_dataset_id: input.preparedDatasetId,
    gene_size: input.geneSize,
    output_dim: input.outputDim,
    use_drug_structure: input.useDrugStructure,
    batch_size: input.batchSize,
    lr: input.lr
  });
  return payload.job;
}

export async function getJob(jobId: string): Promise<JobRecord> {
  const payload = await requestJson<{ job: JobRecord }>(`/api/jobs/${jobId}`);
  return payload.job;
}

export async function cancelJob(jobId: string): Promise<JobRecord> {
  const payload = await requestJson<{ job: JobRecord }>(
    `/api/jobs/${jobId}/cancel`,
    "POST"
  );
  return payload.job;
}

export async function getJobLog(jobId: string): Promise<string> {
  const payload = await requestJson<{ log: string }>(`/api/jobs/${jobId}/log`);
  return payload.log;
}

export async function getModel(modelId: string): Promise<ModelRecord> {
  const payload = await requestJson<{ model: ModelRecord }>(
    `/api/results/models/${modelId}`
  );
  return payload.model;
}

export async function getResultByJob(jobId: string): Promise<ResultRecord> {
  const payload = await requestJson<{ result: ResultRecord }>(`/api/results/job/${jobId}`);
  return payload.result;
}

export type CondaEnvsResponse = {
  conda_bat_candidates: string[];
  conda_envs: string[];
};

export async function getCondaEnvs(condaBat?: string): Promise<CondaEnvsResponse> {
  const path =
    condaBat != null && condaBat !== ""
      ? `/api/runtime/conda-envs?conda_bat=${encodeURIComponent(condaBat)}`
      : "/api/runtime/conda-envs";
  return requestJson<CondaEnvsResponse>(path);
}
