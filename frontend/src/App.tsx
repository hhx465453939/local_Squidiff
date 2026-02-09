import { useEffect, useMemo, useState } from "react";
import {
  buildApiUrl,
  createTrainJob,
  DatasetRecord,
  getHealth,
  getJob,
  getJobLog,
  getModel,
  getResultByJob,
  inspectSeurat,
  JobRecord,
  ModelRecord,
  prepareTraining,
  PrepareTrainingResult,
  ResultRecord,
  SeuratInspectReport,
  uploadDataset,
  validateDataset,
  ValidationReport
} from "./services/api";

function formatTime(value?: string | null): string {
  if (!value) {
    return "-";
  }
  return new Date(value).toLocaleString();
}

function parseSelectedClusters(input: string): string[] {
  return input
    .split(/[,\n;，；]/g)
    .map((item) => item.trim())
    .filter((item) => item.length > 0);
}

export function App() {
  const [health, setHealth] = useState("checking");
  const [globalError, setGlobalError] = useState<string | null>(null);

  const [datasetFile, setDatasetFile] = useState<File | null>(null);
  const [smilesFile, setSmilesFile] = useState<File | null>(null);
  const [datasetName, setDatasetName] = useState("");
  const [useDrugStructure, setUseDrugStructure] = useState(false);

  const [dataset, setDataset] = useState<DatasetRecord | null>(null);
  const [validation, setValidation] = useState<ValidationReport | null>(null);
  const [seuratInspect, setSeuratInspect] = useState<SeuratInspectReport | null>(null);
  const [prepareSummary, setPrepareSummary] = useState<PrepareTrainingResult | null>(null);
  const [groupColumn, setGroupColumn] = useState("");
  const [clusterColumn, setClusterColumn] = useState("");
  const [selectedClustersText, setSelectedClustersText] = useState("");
  const [prepareSeed, setPrepareSeed] = useState(42);
  const [rExecMode, setRExecMode] = useState<"direct" | "cmd_conda">("direct");
  const [rCondaEnv, setRCondaEnv] = useState("");
  const [rCondaBat, setRCondaBat] = useState("conda.bat");
  const [rscriptBin, setRscriptBin] = useState("Rscript");

  const [geneSize, setGeneSize] = useState<number>(100);
  const [outputDim, setOutputDim] = useState<number>(100);
  const [batchSize, setBatchSize] = useState<number>(64);
  const [lr, setLr] = useState<number>(1e-4);

  const [job, setJob] = useState<JobRecord | null>(null);
  const [jobLog, setJobLog] = useState<string>("");
  const [model, setModel] = useState<ModelRecord | null>(null);
  const [result, setResult] = useState<ResultRecord | null>(null);

  const [busyStep, setBusyStep] = useState<string | null>(null);

  useEffect(() => {
    getHealth()
      .then((res) => setHealth(res.status))
      .catch((err: unknown) => {
        setHealth("down");
        setGlobalError(err instanceof Error ? err.message : String(err));
      });
  }, []);

  const canValidate = useMemo(() => dataset !== null, [dataset]);
  const parsedSelectedClusters = useMemo(
    () => parseSelectedClusters(selectedClustersText),
    [selectedClustersText]
  );
  const canPrepare = useMemo(
    () =>
      Boolean(dataset?.path_h5ad) &&
      groupColumn.trim().length > 0 &&
      clusterColumn.trim().length > 0 &&
      parsedSelectedClusters.length > 0,
    [clusterColumn, dataset?.path_h5ad, groupColumn, parsedSelectedClusters.length]
  );
  const canTrain = useMemo(
    () => dataset !== null && validation?.valid === true,
    [dataset, validation]
  );

  useEffect(() => {
    if (job === null) {
      return;
    }
    if (job.status !== "queued" && job.status !== "running") {
      return;
    }
    const timer = window.setTimeout(async () => {
      try {
        const latest = await getJob(job.id);
        setJob(latest);
      } catch (err: unknown) {
        setGlobalError(err instanceof Error ? err.message : String(err));
      }
    }, 2000);
    return () => window.clearTimeout(timer);
  }, [job]);

  useEffect(() => {
    if (job?.status !== "success" && job?.status !== "failed") {
      return;
    }
    getJobLog(job.id)
      .then(setJobLog)
      .catch(() => setJobLog("No log available"));

    if (job.model_id) {
      getModel(job.model_id)
        .then(setModel)
        .catch(() => setModel(null));
    }
    if (job.result_id) {
      getResultByJob(job.id)
        .then(setResult)
        .catch(() => setResult(null));
    }
  }, [job]);

  async function onUpload() {
    if (!datasetFile) {
      setGlobalError("请先选择数据文件");
      return;
    }
    setBusyStep("upload");
    setGlobalError(null);
    try {
      const uploaded = await uploadDataset({
        datasetFile,
        datasetName: datasetName || undefined,
        smilesFile
      });
      setDataset(uploaded);
      setValidation(null);
      setSeuratInspect(null);
      setPrepareSummary(null);
      setGroupColumn("");
      setClusterColumn("");
      setSelectedClustersText("");
      setPrepareSeed(42);
      setJob(null);
      setModel(null);
      setResult(null);
      setJobLog("");
    } catch (err: unknown) {
      setGlobalError(err instanceof Error ? err.message : String(err));
    } finally {
      setBusyStep(null);
    }
  }

  async function onValidate() {
    if (!dataset) {
      setGlobalError("请先上传数据");
      return;
    }
    if (rExecMode === "cmd_conda" && !rCondaEnv.trim()) {
      setGlobalError("请选择 cmd_conda 时必须填写 Conda R 环境名");
      return;
    }
    setBusyStep("validate");
    setGlobalError(null);
    try {
      const payload = await validateDataset({
        datasetId: dataset.id,
        useDrugStructure,
        rExecMode,
        rCondaEnv: rCondaEnv || undefined,
        rCondaBat: rCondaBat || undefined,
        rscriptBin: rscriptBin || undefined
      });
      setDataset(payload.dataset);
      setValidation(payload.validation);
      if (payload.validation.summary?.n_genes) {
        setGeneSize(payload.validation.summary.n_genes);
        setOutputDim(payload.validation.summary.n_genes);
      }
    } catch (err: unknown) {
      setGlobalError(err instanceof Error ? err.message : String(err));
    } finally {
      setBusyStep(null);
    }
  }

  async function onInspectSeurat() {
    if (!dataset) {
      setGlobalError("请先上传数据");
      return;
    }

    setBusyStep("inspect");
    setGlobalError(null);
    try {
      const report = await inspectSeurat({
        datasetId: dataset.id,
        umapPreviewLimit: 500
      });
      setSeuratInspect(report);
      if (!groupColumn.trim()) {
        if (report.metadata_columns.includes("Group")) {
          setGroupColumn("Group");
        } else if (report.metadata_columns.length > 0) {
          setGroupColumn(report.metadata_columns[0]);
        }
      }
      if (!clusterColumn.trim()) {
        if (report.metadata_columns.includes("Cluster")) {
          setClusterColumn("Cluster");
        } else if (report.metadata_columns.includes("celltype")) {
          setClusterColumn("celltype");
        } else if (report.metadata_columns.length > 1) {
          setClusterColumn(report.metadata_columns[1]);
        } else if (report.metadata_columns.length > 0) {
          setClusterColumn(report.metadata_columns[0]);
        }
      }
    } catch (err: unknown) {
      setGlobalError(err instanceof Error ? err.message : String(err));
    } finally {
      setBusyStep(null);
    }
  }

  async function onPrepareTraining() {
    if (!dataset) {
      setGlobalError("请先上传并校验数据");
      return;
    }
    if (!canPrepare) {
      setGlobalError("请填写分组字段、类群字段和待筛选 clusters");
      return;
    }
    setBusyStep("prepare");
    setGlobalError(null);
    try {
      const prepared = await prepareTraining({
        datasetId: dataset.id,
        groupColumn: groupColumn.trim(),
        clusterColumn: clusterColumn.trim(),
        selectedClusters: parsedSelectedClusters,
        seed: prepareSeed
      });
      setPrepareSummary(prepared);
      setGeneSize(prepared.n_genes);
      setOutputDim(prepared.n_genes);
    } catch (err: unknown) {
      setGlobalError(err instanceof Error ? err.message : String(err));
    } finally {
      setBusyStep(null);
    }
  }

  async function onSubmitTrain() {
    if (!dataset) {
      setGlobalError("请先上传并校验数据");
      return;
    }
    setBusyStep("train");
    setGlobalError(null);
    try {
      const newJob = await createTrainJob({
        datasetId: dataset.id,
        preparedDatasetId: prepareSummary?.prepared_dataset_id,
        geneSize,
        outputDim,
        useDrugStructure,
        batchSize,
        lr
      });
      setJob(newJob);
      setModel(null);
      setResult(null);
      setJobLog("");
    } catch (err: unknown) {
      setGlobalError(err instanceof Error ? err.message : String(err));
    } finally {
      setBusyStep(null);
    }
  }

  return (
    <main className="app-shell">
      <header className="hero">
        <h1>Squidiff LabFlow MVP</h1>
        <p>上传 → 校验 → Seurat 检查 → 500x500 预处理 → 训练任务 → 结果页</p>
      </header>

      <section className="panel">
        <h2>系统状态</h2>
        <p className="status">backend: {health}</p>
        {globalError ? <p className="error">{globalError}</p> : null}
      </section>

      <section className="panel">
        <h2>1) 上传数据</h2>
        <div className="form-grid">
          <label>
            数据名称
            <input
              value={datasetName}
              onChange={(e) => setDatasetName(e.target.value)}
              placeholder="可选"
            />
          </label>
          <label>
            主数据文件（h5ad/rds/h5seurat）
            <input
              type="file"
              accept=".h5ad,.rds,.h5seurat"
              onChange={(e) => setDatasetFile(e.target.files?.[0] ?? null)}
            />
          </label>
          <label>
            SMILES CSV（可选）
            <input
              type="file"
              accept=".csv,.tsv,.txt"
              onChange={(e) => setSmilesFile(e.target.files?.[0] ?? null)}
            />
          </label>
          <label className="checkbox">
            <input
              type="checkbox"
              checked={useDrugStructure}
              onChange={(e) => setUseDrugStructure(e.target.checked)}
            />
            使用药物结构模式（需要 SMILES + dose 字段）
          </label>
        </div>
        <button disabled={busyStep !== null} onClick={onUpload}>
          {busyStep === "upload" ? "上传中..." : "上传"}
        </button>
        {dataset ? (
          <div className="kv">
            <span>dataset_id</span>
            <code>{dataset.id}</code>
            <span>status</span>
            <code>{dataset.status}</code>
          </div>
        ) : null}
      </section>

      <section className="panel">
        <h2>2) 校验</h2>
        <div className="form-grid compact">
          <label>
            R 执行方式
            <select
              value={rExecMode}
              onChange={(e) => setRExecMode(e.target.value as "direct" | "cmd_conda")}
            >
              <option value="direct">direct (直接调用 Rscript)</option>
              <option value="cmd_conda">cmd_conda (cmd + conda activate)</option>
            </select>
          </label>
          <label>
            Rscript 命令
            <input
              value={rscriptBin}
              onChange={(e) => setRscriptBin(e.target.value)}
              placeholder="Rscript"
            />
          </label>
          <label>
            conda.bat 路径
            <input
              value={rCondaBat}
              onChange={(e) => setRCondaBat(e.target.value)}
              placeholder="conda.bat 或完整路径"
            />
          </label>
          <label>
            Conda R 环境名
            <input
              value={rCondaEnv}
              onChange={(e) => setRCondaEnv(e.target.value)}
              placeholder="例如: r-seurat"
            />
          </label>
        </div>
        <button disabled={!canValidate || busyStep !== null} onClick={onValidate}>
          {busyStep === "validate" ? "校验中..." : "开始校验"}
        </button>
        {validation ? (
          <div className="stack">
            <p>
              校验结果:
              <strong>{validation.valid ? "通过" : "失败"}</strong>
            </p>
            {validation.summary ? (
              <div className="kv">
                <span>cells</span>
                <code>{validation.summary.n_cells}</code>
                <span>genes</span>
                <code>{validation.summary.n_genes}</code>
              </div>
            ) : null}
            {validation.errors.length > 0 ? (
              <ul className="error-list">
                {validation.errors.map((item) => (
                  <li key={item}>{item}</li>
                ))}
              </ul>
            ) : null}
            {validation.warnings.length > 0 ? (
              <ul>
                {validation.warnings.map((item) => (
                  <li key={item}>{item}</li>
                ))}
              </ul>
            ) : null}
          </div>
        ) : null}
      </section>

      <section className="panel">
        <h2>3) Seurat 解析（V2 Phase 1）</h2>
        <p>读取 metadata 字段和 UMAP 预览（rds/h5seurat 需先校验转换为 h5ad）。</p>
        <button
          disabled={!dataset?.path_h5ad || busyStep !== null}
          onClick={onInspectSeurat}
        >
          {busyStep === "inspect" ? "解析中..." : "解析 Seurat"}
        </button>

        {seuratInspect ? (
          <div className="stack">
            <div className="kv">
              <span>n_cells</span>
              <code>{seuratInspect.n_cells}</code>
              <span>n_genes</span>
              <code>{seuratInspect.n_genes}</code>
              <span>has_umap</span>
              <code>{seuratInspect.has_umap ? "yes" : "no"}</code>
            </div>

            <div>
              <h3>metadata 字段</h3>
              {seuratInspect.metadata_columns.length > 0 ? (
                <div className="chip-list">
                  {seuratInspect.metadata_columns.map((column) => (
                    <code key={column}>{column}</code>
                  ))}
                </div>
              ) : (
                <p>未找到 metadata 列。</p>
              )}
            </div>

            <div>
              <h3>UMAP 预览</h3>
              {seuratInspect.umap ? (
                <div className="stack">
                  <div className="kv">
                    <span>key</span>
                    <code>{seuratInspect.umap.key}</code>
                    <span>points</span>
                    <code>
                      {seuratInspect.umap.preview_count}/{seuratInspect.umap.n_points}
                    </code>
                  </div>
                  <pre className="log umap-preview">
                    {JSON.stringify(seuratInspect.umap.preview.slice(0, 30), null, 2)}
                  </pre>
                </div>
              ) : (
                <p>无 UMAP 数据（不阻断后续流程）。</p>
              )}
            </div>

            {seuratInspect.warnings.length > 0 ? (
              <ul>
                {seuratInspect.warnings.map((item) => (
                  <li key={item}>{item}</li>
                ))}
              </ul>
            ) : null}
          </div>
        ) : null}
      </section>

      <section className="panel">
        <h2>4) 500x500 预处理（V2 Phase 2）</h2>
        <div className="form-grid compact">
          <label>
            分组字段（group_column）
            <input
              value={groupColumn}
              onChange={(e) => setGroupColumn(e.target.value)}
              placeholder="例如: sample"
              list="metadata-columns"
            />
          </label>
          <label>
            类群字段（cluster_column）
            <input
              value={clusterColumn}
              onChange={(e) => setClusterColumn(e.target.value)}
              placeholder="例如: celltype"
              list="metadata-columns"
            />
          </label>
          <label>
            随机种子（seed）
            <input
              type="number"
              value={prepareSeed}
              onChange={(e) => setPrepareSeed(Number(e.target.value))}
            />
          </label>
          <label>
            待筛选 clusters（逗号分隔）
            <input
              value={selectedClustersText}
              onChange={(e) => setSelectedClustersText(e.target.value)}
              placeholder="例如: T,B,NK"
            />
          </label>
        </div>
        <datalist id="metadata-columns">
          {seuratInspect?.metadata_columns.map((column) => (
            <option key={column} value={column} />
          ))}
        </datalist>
        <button disabled={!canPrepare || busyStep !== null} onClick={onPrepareTraining}>
          {busyStep === "prepare" ? "预处理中..." : "准备训练数据（500x500）"}
        </button>
        {prepareSummary ? (
          <div className="stack">
            <p className="status">预处理完成，可直接进入训练。</p>
            <div className="kv">
              <span>prepared_dataset_id</span>
              <code>{prepareSummary.prepared_dataset_id}</code>
              <span>n_cells</span>
              <code>{prepareSummary.n_cells}</code>
              <span>n_genes</span>
              <code>{prepareSummary.n_genes}</code>
              <span>sampling_method</span>
              <code>{prepareSummary.sampling_report.mode}</code>
              <span>gene_method</span>
              <code>{prepareSummary.gene_report.method}</code>
            </div>
          </div>
        ) : (
          <p>未执行预处理时，训练接口会自动尝试使用该数据集最近一次 prepared dataset。</p>
        )}
      </section>

      <section className="panel">
        <h2>5) 提交训练任务</h2>
        {prepareSummary ? (
          <p className="status">
            当前训练默认使用 prepared_dataset_id: {prepareSummary.prepared_dataset_id}
          </p>
        ) : (
          <p>当前训练将使用 dataset_id，若存在历史 prepared 数据集将由后端自动选择最新版本。</p>
        )}
        <div className="form-grid compact">
          <label>
            gene_size
            <input
              type="number"
              value={geneSize}
              onChange={(e) => setGeneSize(Number(e.target.value))}
            />
          </label>
          <label>
            output_dim
            <input
              type="number"
              value={outputDim}
              onChange={(e) => setOutputDim(Number(e.target.value))}
            />
          </label>
          <label>
            batch_size
            <input
              type="number"
              value={batchSize}
              onChange={(e) => setBatchSize(Number(e.target.value))}
            />
          </label>
          <label>
            lr
            <input
              type="number"
              step="0.00001"
              value={lr}
              onChange={(e) => setLr(Number(e.target.value))}
            />
          </label>
        </div>
        <button disabled={!canTrain || busyStep !== null} onClick={onSubmitTrain}>
          {busyStep === "train" ? "提交中..." : "开始训练"}
        </button>
      </section>

      <section className="panel">
        <h2>6) 任务轮询状态</h2>
        {job ? (
          <div className="stack">
            <div className="kv">
              <span>job_id</span>
              <code>{job.id}</code>
              <span>status</span>
              <code>{job.status}</code>
              <span>source_dataset_id</span>
              <code>{job.source_dataset_id ?? "-"}</code>
              <span>train_dataset_id</span>
              <code>{job.dataset_id}</code>
              <span>prepared_dataset_id</span>
              <code>{job.prepared_dataset_id ?? "-"}</code>
              <span>started</span>
              <code>{formatTime(job.started_at)}</code>
              <span>ended</span>
              <code>{formatTime(job.ended_at)}</code>
            </div>
            {job.error_msg ? <pre className="log">{job.error_msg}</pre> : null}
          </div>
        ) : (
          <p>尚未创建任务</p>
        )}
      </section>

      <section className="panel">
        <h2>7) 结果页</h2>
        {job?.status === "success" ? (
          <div className="stack">
            <p>训练已完成。</p>
            {model ? (
              <div className="kv">
                <span>model_id</span>
                <code>{model.id}</code>
                <span>checkpoint</span>
                <code className="path">{model.path_ckpt}</code>
                <span>gene_size</span>
                <code>{model.gene_size}</code>
                <span>output_dim</span>
                <code>{model.output_dim}</code>
              </div>
            ) : (
              <p>模型记录加载中或未生成。</p>
            )}

            {result ? (
              <div className="stack">
                <h3>预测结果</h3>
                <div className="asset-grid">
                  {result.summary.assets.map((asset) => (
                    <figure key={asset.name} className="asset-card">
                      <img src={buildApiUrl(asset.url)} alt={asset.type} />
                      <figcaption>{asset.type}</figcaption>
                    </figure>
                  ))}
                </div>
              </div>
            ) : null}

            {jobLog ? (
              <details>
                <summary>查看日志</summary>
                <pre className="log">{jobLog}</pre>
              </details>
            ) : null}
          </div>
        ) : (
          <p>训练完成后会显示模型与结果信息。</p>
        )}
      </section>
    </main>
  );
}
