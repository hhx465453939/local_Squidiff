import { useCallback, useEffect, useMemo, useState } from "react";
import { ParamTooltip, PARAM_TOOLTIPS } from "./components/ParamTooltip";
import {
  AuthUser,
  buildApiUrl,
  cancelJob,
  createTrainJob,
  DatasetRecord,
  deleteJob,
  deleteModel,
  deleteResult,
  flushJobs,
  getCurrentUser,
  getCondaEnvs,
  getGpuStats,
  getHealth,
  getJob,
  getJobs,
  getJobLog,
  getSchedulerPref,
  GpuStatsResponse,
  listModels,
  listResults,
  loginUser,
  logoutUser,
  getModel,
  getResultByJob,
  inspectSeurat,
  JobRecord,
  ModelRecord,
  prepareTraining,
  PrepareTrainingResult,
  registerUser,
  ResultRecord,
  SeuratInspectReport,
  setAuthToken,
  updateSchedulerPref,
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

function clampPercent(value: number): number {
  if (!Number.isFinite(value)) {
    return 0;
  }
  return Math.max(0, Math.min(100, value));
}

const ACTIVE_JOB_STATUS = new Set<JobRecord["status"]>(["queued", "running"]);

function getJobSortTime(job: JobRecord): number {
  const raw = job.updated_at ?? job.created_at;
  const parsed = Date.parse(raw);
  return Number.isNaN(parsed) ? 0 : parsed;
}

function sortJobsByRecent(items: JobRecord[]): JobRecord[] {
  return [...items].sort((a, b) => getJobSortTime(b) - getJobSortTime(a));
}

function mergeJobsById(current: JobRecord[], incoming: JobRecord[]): JobRecord[] {
  const merged = new Map<string, JobRecord>();
  for (const job of current) {
    merged.set(job.id, job);
  }
  for (const job of incoming) {
    merged.set(job.id, job);
  }
  return sortJobsByRecent(Array.from(merged.values()));
}

function sortByUpdated<T extends { updated_at: string; created_at: string }>(
  items: T[]
): T[] {
  return [...items].sort((a, b) => {
    const at = Date.parse(a.updated_at ?? a.created_at);
    const bt = Date.parse(b.updated_at ?? b.created_at);
    const aa = Number.isNaN(at) ? 0 : at;
    const bb = Number.isNaN(bt) ? 0 : bt;
    return bb - aa;
  });
}

export function App() {
  const [health, setHealth] = useState("checking");
  const [globalError, setGlobalError] = useState<string | null>(null);
  const [authUser, setAuthUser] = useState<AuthUser | null>(null);
  const [authMode, setAuthMode] = useState<"login" | "register">("login");
  const [authUsername, setAuthUsername] = useState("");
  const [authPassword, setAuthPassword] = useState("");
  const [authBusy, setAuthBusy] = useState(false);
  const [authError, setAuthError] = useState<string | null>(null);

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
  const [validateError, setValidateError] = useState<string | null>(null);
  const [condaBatCandidates, setCondaBatCandidates] = useState<string[]>([]);
  const [condaEnvsList, setCondaEnvsList] = useState<string[]>([]);

  const [geneSize, setGeneSize] = useState<number>(100);
  const [outputDim, setOutputDim] = useState<number>(100);
  const [batchSize, setBatchSize] = useState<number>(64);
  const [lr, setLr] = useState<number>(1e-4);

  const [job, setJob] = useState<JobRecord | null>(null);
  const [jobHistory, setJobHistory] = useState<JobRecord[]>([]);
  const [modelHistory, setModelHistory] = useState<ModelRecord[]>([]);
  const [resultHistory, setResultHistory] = useState<ResultRecord[]>([]);
  const [jobHistoryLoading, setJobHistoryLoading] = useState(false);
  const [jobLog, setJobLog] = useState<string>("");
  const [liveJobLog, setLiveJobLog] = useState<string>("");
  const [gpuStats, setGpuStats] = useState<GpuStatsResponse | null>(null);
  const [gpuStatsError, setGpuStatsError] = useState<string | null>(null);
  const [model, setModel] = useState<ModelRecord | null>(null);
  const [result, setResult] = useState<ResultRecord | null>(null);
  const [schedulerMode, setSchedulerMode] = useState<"serial" | "parallel">("parallel");
  const [schedulerLimit, setSchedulerLimit] = useState(3);
  const [schedulerBusy, setSchedulerBusy] = useState(false);
  const [schedulerError, setSchedulerError] = useState<string | null>(null);

  const [busyStep, setBusyStep] = useState<string | null>(null);
  const [cancelBusy, setCancelBusy] = useState(false);
  const [hasEntered, setHasEntered] = useState(false);
  const [selectedMetadataColumn, setSelectedMetadataColumn] = useState<string | null>(null);
  const jobSessionKey = authUser ? `labflow_last_job_${authUser.id}` : null;

  function resetCurrentJobView(): void {
    setJob(null);
    setModel(null);
    setResult(null);
    setJobLog("");
    setLiveJobLog("");
    setGpuStats(null);
    setGpuStatsError(null);
  }

  function onSelectJob(jobItem: JobRecord): void {
    setJob(jobItem);
    setModel(null);
    setResult(null);
    setJobLog("");
    setLiveJobLog("");
    setGlobalError(null);
  }

  function onStartNewTask(): void {
    setDatasetFile(null);
    setSmilesFile(null);
    setDatasetName("");
    setUseDrugStructure(false);
    setDataset(null);
    setValidation(null);
    setSeuratInspect(null);
    setPrepareSummary(null);
    setGroupColumn("");
    setClusterColumn("");
    setSelectedClustersText("");
    setSelectedMetadataColumn(null);
    setPrepareSeed(42);
    setGeneSize(100);
    setOutputDim(100);
    setBatchSize(64);
    setLr(1e-4);
    resetCurrentJobView();
    if (jobSessionKey && typeof window !== "undefined") {
      window.localStorage.removeItem(jobSessionKey);
    }
  }

  useEffect(() => {
    getHealth()
      .then((res) => setHealth(res.status))
      .catch((err: unknown) => {
        setHealth("down");
        setGlobalError(err instanceof Error ? err.message : String(err));
      });
  }, []);

  useEffect(() => {
    getCurrentUser()
      .then((user) => setAuthUser(user))
      .catch(() => setAuthUser(null));
  }, []);

  useEffect(() => {
    if (!authUser) {
      return;
    }
    getCondaEnvs()
      .then((res) => {
        const uniqueCandidates = Array.from(new Set(res.conda_bat_candidates));
        setCondaBatCandidates(uniqueCandidates);
        // 如果后端返回了环境列表，直接使用
        if (res.conda_envs.length > 0) {
          const uniqueEnvs = Array.from(new Set(res.conda_envs));
          setCondaEnvsList(uniqueEnvs);
          const rEnv = uniqueEnvs.find((e) => /^r[-.]?\d|r-seurat/i.test(e));
          if (rEnv) {
            setRCondaEnv((prev) => (prev === "" ? rEnv : prev));
          }
        }
        // 如果有 conda.bat 候选但环境列表为空，用第一个 candidate 重新获取
        if (uniqueCandidates.length > 0) {
          const firstBat = uniqueCandidates[0];
          setRCondaBat((prev) =>
            prev === "conda.bat" || !uniqueCandidates.includes(prev)
              ? firstBat
              : prev
          );
          // 如果环境列表为空，用第一个 candidate 重新获取
          if (res.conda_envs.length === 0) {
            getCondaEnvs(firstBat)
              .then((res2) => {
                const uniqueEnvs = Array.from(new Set(res2.conda_envs));
                setCondaEnvsList(uniqueEnvs);
                const rEnv = uniqueEnvs.find((e) => /^r[-.]?\d|r-seurat/i.test(e));
                if (rEnv) {
                  setRCondaEnv((prev) => (prev === "" ? rEnv : prev));
                }
              })
              .catch(() => {});
          }
        }
      })
      .catch(() => {});
  }, [authUser]);

  const loadJobHistory = useCallback(
    async (options?: { selectPreferred?: boolean; showLoading?: boolean }) => {
      const selectPreferred = options?.selectPreferred ?? false;
      const showLoading = options?.showLoading ?? false;

      if (!authUser) {
        setJobHistory([]);
        return;
      }

      if (showLoading) {
        setJobHistoryLoading(true);
      }
      try {
        const [jobs, models, results] = await Promise.all([
          getJobs(),
          listModels(),
          listResults()
        ]);
        const items = sortJobsByRecent(jobs);
        setJobHistory(items);
        setModelHistory(sortByUpdated(models));
        setResultHistory(sortByUpdated(results));

        setJob((current) => {
          const matchedCurrent = current
            ? items.find((item) => item.id === current.id) ?? null
            : null;
          if (matchedCurrent) {
            return matchedCurrent;
          }
          if (!selectPreferred) {
            return null;
          }

          const storedJobId =
            jobSessionKey && typeof window !== "undefined"
              ? window.localStorage.getItem(jobSessionKey)
              : null;
          const preferred =
            (storedJobId ? items.find((item) => item.id === storedJobId) : null) ??
            items.find((item) => ACTIVE_JOB_STATUS.has(item.status)) ??
            items[0] ??
            null;
          if (preferred) {
            setHasEntered(true);
          }
          return preferred;
        });
      } catch (err: unknown) {
        setGlobalError(err instanceof Error ? err.message : String(err));
      } finally {
        if (showLoading) {
          setJobHistoryLoading(false);
        }
      }
    },
    [authUser, jobSessionKey]
  );
  const selectedJobId = job?.id ?? null;
  const recentJobs = useMemo(() => jobHistory.slice(0, 20), [jobHistory]);
  const recentModels = useMemo(() => modelHistory.slice(0, 20), [modelHistory]);
  const recentResults = useMemo(() => resultHistory.slice(0, 20), [resultHistory]);

  useEffect(() => {
    if (!authUser) {
      setJobHistory([]);
      setModelHistory([]);
      setResultHistory([]);
      setJob(null);
      setModel(null);
      setResult(null);
      setJobLog("");
      setLiveJobLog("");
      setGpuStats(null);
      setGpuStatsError(null);
      setSchedulerMode("parallel");
      setSchedulerLimit(3);
      setSchedulerError(null);
      return;
    }
    loadJobHistory({ selectPreferred: true, showLoading: true }).catch(() => {});
  }, [authUser, loadJobHistory]);

  useEffect(() => {
    if (!authUser || !hasEntered) {
      return;
    }
    const tick = () => {
      loadJobHistory({ selectPreferred: false, showLoading: false }).catch(() => {});
    };
    tick();
    const interval = window.setInterval(tick, 5000);
    return () => window.clearInterval(interval);
  }, [authUser, hasEntered, job?.id, loadJobHistory]);

  useEffect(() => {
    if (!authUser) {
      return;
    }
    getSchedulerPref()
      .then((pref) => {
        setSchedulerMode(pref.mode);
        setSchedulerLimit(pref.max_concurrent_jobs);
        setSchedulerError(null);
      })
      .catch((err: unknown) => {
        setSchedulerError(err instanceof Error ? err.message : String(err));
      });
  }, [authUser]);

  useEffect(() => {
    if (!jobSessionKey || typeof window === "undefined") {
      return;
    }
    if (!job) {
      window.localStorage.removeItem(jobSessionKey);
      return;
    }
    window.localStorage.setItem(jobSessionKey, job.id);
  }, [jobSessionKey, job]);

  const canValidate = useMemo(() => dataset !== null, [dataset]);
  const activeJobs = useMemo(
    () => jobHistory.filter((item) => ACTIVE_JOB_STATUS.has(item.status)),
    [jobHistory]
  );
  const myActiveJobsCount = useMemo(() => {
    if (!authUser) {
      return 0;
    }
    return activeJobs.filter((item) => item.owner_user_id === authUser.id).length;
  }, [activeJobs, authUser]);
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
    // 训练按钮放宽：校验通过 or 500×500 预处理完成 均允许开始训练
    () => dataset !== null && (validation?.valid === true || prepareSummary !== null),
    [dataset, prepareSummary, validation]
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
        setJobHistory((prev) => mergeJobsById(prev, [latest]));
      } catch (err: unknown) {
        setGlobalError(err instanceof Error ? err.message : String(err));
      }
    }, 2000);
    return () => window.clearTimeout(timer);
  }, [job]);

  // 训练/任务运行中时轮询运行日志，实现「小电视」实时反馈
  useEffect(() => {
    if (job == null || (job.status !== "queued" && job.status !== "running")) {
      return;
    }
    const jobId = job.id;
    const tick = () => {
      getJobLog(jobId)
        .then(setLiveJobLog)
        .catch(() => {});
    };
    tick();
    const interval = window.setInterval(tick, 2500);
    return () => window.clearInterval(interval);
  }, [job]);

  useEffect(() => {
    if (job?.status !== "running") {
      setGpuStats(null);
      setGpuStatsError(null);
      return;
    }
    const tick = () => {
      getGpuStats()
        .then((stats) => {
          setGpuStats(stats);
          setGpuStatsError(null);
        })
        .catch((err: unknown) => {
          setGpuStatsError(err instanceof Error ? err.message : String(err));
        });
    };
    tick();
    const interval = window.setInterval(tick, 3000);
    return () => window.clearInterval(interval);
  }, [job?.id, job?.status]);

  useEffect(() => {
    if (
      job?.status !== "success" &&
      job?.status !== "failed" &&
      job?.status !== "canceled"
    ) {
      return;
    }
    getJobLog(job.id)
      .then((log) => {
        setJobLog(log);
        setLiveJobLog(log);
      })
      .catch(() => {
        setJobLog("No log available");
        setLiveJobLog("");
      });

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
      setLiveJobLog("");
    } catch (err: unknown) {
      setGlobalError(err instanceof Error ? err.message : String(err));
    } finally {
      setBusyStep(null);
    }
  }

  async function onValidate() {
    if (!dataset) {
      setValidateError("请先上传数据");
      return;
    }
    if (rExecMode === "cmd_conda" && !rCondaEnv.trim()) {
      setValidateError("请选择 cmd_conda 时必须填写 Conda R 环境名");
      return;
    }
    setBusyStep("validate");
    setValidateError(null);
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
      setValidateError(null);
      if (payload.validation.summary?.n_genes) {
        setGeneSize(payload.validation.summary.n_genes);
        setOutputDim(payload.validation.summary.n_genes);
      }
    } catch (err: unknown) {
      const raw = err instanceof Error ? err.message : String(err);
      let msg = raw;
      try {
        const parsed = JSON.parse(raw) as { detail?: string };
        if (typeof parsed.detail === "string") {
          msg = parsed.detail;
        }
      } catch {
        // not JSON, keep raw
      }
      setValidateError(msg);
    } finally {
      setBusyStep(null);
    }
  }

  const validateRecommendation =
    validateError && /Rscript|cmd_conda|Conda|conda\.bat/i.test(validateError)
      ? "推荐：R 执行方式选 cmd_conda，conda.bat 填完整路径（如 F:\\software\\Miniconda3\\condabin\\conda.bat），Conda R 环境名填 r-4.3"
      : null;

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
      setJobHistory((prev) => mergeJobsById(prev, [newJob]));
      setModel(null);
      setResult(null);
      setJobLog("");
      setLiveJobLog("");
    } catch (err: unknown) {
      setGlobalError(err instanceof Error ? err.message : String(err));
    } finally {
      setBusyStep(null);
    }
  }

  async function onCancelTrain() {
    if (!job || (job.status !== "queued" && job.status !== "running")) {
      return;
    }
    setCancelBusy(true);
    setGlobalError(null);
    try {
      const updated = await cancelJob(job.id);
      setJob(updated);
      setJobHistory((prev) => mergeJobsById(prev, [updated]));
      const log = await getJobLog(job.id);
      setLiveJobLog(log);
      setJobLog(log);
    } catch (err: unknown) {
      setGlobalError(err instanceof Error ? err.message : String(err));
    } finally {
      setCancelBusy(false);
    }
  }

  async function onAuthSubmit() {
    const username = authUsername.trim();
    const password = authPassword;
    if (!username || !password) {
      setAuthError("请输入用户名和密码");
      return;
    }
    setAuthBusy(true);
    setAuthError(null);
    try {
      const payload =
        authMode === "register"
          ? await registerUser({ username, password })
          : await loginUser({ username, password });
      setAuthToken(payload.access_token);
      setAuthUser(payload.user);
      setAuthPassword("");
      setHasEntered(true);
    } catch (err: unknown) {
      const raw = err instanceof Error ? err.message : String(err);
      let msg = raw;
      try {
        const parsed = JSON.parse(raw) as { detail?: string };
        if (typeof parsed.detail === "string") {
          msg = parsed.detail;
        }
      } catch {
        // keep raw
      }
      setAuthError(msg);
    } finally {
      setAuthBusy(false);
    }
  }

  async function onAuthLogout() {
    setAuthBusy(true);
    setAuthError(null);
    try {
      await logoutUser();
    } catch {
      // ignore remote logout failure, still clear local token
    } finally {
      setAuthToken(null);
      setAuthUser(null);
      setHasEntered(false);
      setJobHistory([]);
      setModelHistory([]);
      setResultHistory([]);
      setJobHistoryLoading(false);
      setAuthBusy(false);
    }
  }

  async function onRefreshJobs() {
    await loadJobHistory({ selectPreferred: false, showLoading: true });
  }

  async function onToggleSchedulerMode(mode: "serial" | "parallel") {
    setSchedulerBusy(true);
    setSchedulerError(null);
    try {
      const pref = await updateSchedulerPref(mode);
      setSchedulerMode(pref.mode);
      setSchedulerLimit(pref.max_concurrent_jobs);
    } catch (err: unknown) {
      setSchedulerError(err instanceof Error ? err.message : String(err));
    } finally {
      setSchedulerBusy(false);
    }
  }

  async function onFlushActiveJobs() {
    const ok = window.confirm(
      "One-click clear all queued/running tasks and related artifacts? This cannot be undone."
    );
    if (!ok) {
      return;
    }
    try {
      const summary = await flushJobs({
        scope: "active",
        purgeArtifacts: true,
        force: true
      });
      await loadJobHistory({ selectPreferred: false, showLoading: true });
      if (summary.skipped_running_jobs.length > 0) {
        window.alert(
          `Cleared ${summary.deleted_jobs} tasks, but ${summary.skipped_running_jobs.length} are still running.`
        );
      } else {
        window.alert(
          `Cleared ${summary.deleted_jobs} tasks (${summary.removed_models} models, ${summary.removed_results} reports).`
        );
      }
    } catch (err: unknown) {
      setGlobalError(err instanceof Error ? err.message : String(err));
    }
  }

  async function openJobById(jobId: string) {
    try {
      const existing = jobHistory.find((item) => item.id === jobId);
      if (existing) {
        onSelectJob(existing);
        return;
      }
      const loaded = await getJob(jobId);
      onSelectJob(loaded);
      setJobHistory((prev) => mergeJobsById(prev, [loaded]));
    } catch (err: unknown) {
      setGlobalError(err instanceof Error ? err.message : String(err));
    }
  }

  async function onDeleteJob(jobId: string) {
    const ok = window.confirm(
      `Delete task ${jobId}? This will remove task metadata and related artifacts.`
    );
    if (!ok) {
      return;
    }
    try {
      await deleteJob(jobId, true);
      if (job?.id === jobId) {
        resetCurrentJobView();
      }
      await loadJobHistory({ selectPreferred: false, showLoading: true });
    } catch (err: unknown) {
      setGlobalError(err instanceof Error ? err.message : String(err));
    }
  }

  async function onDeleteModel(modelId: string) {
    const ok = window.confirm(
      `Delete model ${modelId}? This will remove the model record and checkpoint file.`
    );
    if (!ok) {
      return;
    }
    try {
      await deleteModel(modelId, true);
      await loadJobHistory({ selectPreferred: false, showLoading: true });
    } catch (err: unknown) {
      setGlobalError(err instanceof Error ? err.message : String(err));
    }
  }

  async function onDeleteResult(resultId: string) {
    const ok = window.confirm(
      `Delete result ${resultId}? This will remove result metadata and report files.`
    );
    if (!ok) {
      return;
    }
    try {
      await deleteResult(resultId, true);
      await loadJobHistory({ selectPreferred: false, showLoading: true });
    } catch (err: unknown) {
      setGlobalError(err instanceof Error ? err.message : String(err));
    }
  }

  const taskLabelMap: Record<string, string> = {
    upload: "上传中…",
    validate: "校验中…",
    inspect: "解析 Seurat 中…",
    prepare: "500×500 预处理中…",
    train: "提交训练中…"
  };
  const isJobActive =
    job != null && (job.status === "queued" || job.status === "running");
  const canCancelTrain = job?.status === "running" && !cancelBusy;
  const showTaskFeedback = busyStep !== null || isJobActive;
  const taskLabel = busyStep
    ? taskLabelMap[busyStep] ?? "处理中…"
    : isJobActive
      ? "训练运行中…"
      : "";

  return (
    <main className={`app-shell ${hasEntered ? "app-shell-entered" : "app-shell-landing"}`}>
      {!hasEntered ? (
        <section className="landing-hero panel">
          <p className="landing-kicker">LabFlow Scientific Workspace</p>
          <h1>Squidiff Research Console</h1>
          <p className="landing-subtitle">
            从单细胞数据上传到训练结果可视化的一体化实验流程。
            点击入口，立即开始本轮分析。
          </p>
          <div className="landing-badges" aria-label="核心能力">
            <span>Seurat Inspect</span>
            <span>500x500 Prepare</span>
            <span>GPU Training</span>
            <span>Live Logs</span>
          </div>
          <div className="landing-grid">
            <article className="landing-card">
              <h2>数据质检优先</h2>
              <p>先校验、后训练，异常在步骤内就近提示，减少反复重跑。</p>
            </article>
            <article className="landing-card">
              <h2>实验过程可观测</h2>
              <p>任务状态、运行日志、模型与结果资产统一归档追踪。</p>
            </article>
            <article className="landing-card">
              <h2>科研节奏友好</h2>
              <p>支持 prepared dataset 复用，快速迭代参数与训练批次。</p>
            </article>
          </div>
          <div className="landing-actions">
            <button
              type="button"
              className="enter-analysis-btn"
              onClick={() => setHasEntered(true)}
              disabled={!authUser}
            >
              开始分析
            </button>
            <a
              className="guide-link-btn"
              href={buildApiUrl("/api/auth/user-guide")}
              target="_blank"
              rel="noreferrer"
            >
              用户说明书
            </a>
            <p className="landing-health">backend: {health}</p>
          </div>
          <div className="auth-panel">
            <p className="auth-title">账户入口</p>
            <div className="auth-tabs" role="tablist" aria-label="登录模式">
              <button
                type="button"
                className={authMode === "login" ? "auth-tab active" : "auth-tab"}
                onClick={() => {
                  setAuthMode("login");
                  setAuthError(null);
                }}
              >
                登录
              </button>
              <button
                type="button"
                className={authMode === "register" ? "auth-tab active" : "auth-tab"}
                onClick={() => {
                  setAuthMode("register");
                  setAuthError(null);
                }}
              >
                注册
              </button>
            </div>
            <div className="auth-fields">
              <label>
                用户名
                <input
                  value={authUsername}
                  onChange={(e) => setAuthUsername(e.target.value)}
                  placeholder="例如: labflow_user"
                />
              </label>
              <label>
                密码
                <input
                  type="password"
                  value={authPassword}
                  onChange={(e) => setAuthPassword(e.target.value)}
                  placeholder="至少 8 位"
                />
              </label>
            </div>
            <div className="auth-actions">
              <button type="button" onClick={onAuthSubmit} disabled={authBusy}>
                {authBusy ? "提交中..." : authMode === "register" ? "注册并登录" : "登录"}
              </button>
              {authUser ? (
                <button
                  type="button"
                  className="auth-ghost-btn"
                  onClick={onAuthLogout}
                  disabled={authBusy}
                >
                  退出登录
                </button>
              ) : null}
            </div>
            <p className="auth-state">
              {authUser ? `当前用户：${authUser.username}` : "未登录，开始分析按钮将保持禁用"}
            </p>
            {authError ? <p className="error">{authError}</p> : null}
          </div>
          {globalError ? <p className="error">{globalError}</p> : null}
        </section>
      ) : (
        <div className="workspace-layout">
      <aside className="panel task-sidebar">
        <h2>Task Center</h2>
        <p className="task-sidebar-meta">
          {jobHistoryLoading
            ? "syncing..."
            : `active ${activeJobs.length} / total ${jobHistory.length}`}
        </p>
        <div className="task-sidebar-actions">
          <button type="button" onClick={onStartNewTask} disabled={busyStep !== null}>
            New task
          </button>
          <button type="button" className="danger-btn" onClick={onFlushActiveJobs}>
            One-click clear
          </button>
          <button type="button" className="auth-ghost-btn" onClick={onRefreshJobs}>
            Refresh
          </button>
          <button type="button" className="auth-ghost-btn" onClick={onAuthLogout}>
            Logout
          </button>
        </div>

        <section className="sidebar-block" aria-label="tasks">
          <h3>Jobs</h3>
          <p className="hint">
            My scheduler: {schedulerMode} ({myActiveJobsCount}/{schedulerLimit})
          </p>
          <div className="auth-tabs" role="tablist" aria-label="scheduler mode">
            <button
              type="button"
              className={schedulerMode === "serial" ? "auth-tab active" : "auth-tab"}
              onClick={() => onToggleSchedulerMode("serial")}
              disabled={schedulerBusy}
            >
              Serial (1)
            </button>
            <button
              type="button"
              className={schedulerMode === "parallel" ? "auth-tab active" : "auth-tab"}
              onClick={() => onToggleSchedulerMode("parallel")}
              disabled={schedulerBusy}
            >
              Parallel (3)
            </button>
          </div>
          {schedulerError ? <p className="error">{schedulerError}</p> : null}
          <div className="task-list">
            {recentJobs.length > 0 ? (
              recentJobs.map((item) => (
                <div
                  key={item.id}
                  className={`task-item ${selectedJobId === item.id ? "active" : ""}`}
                >
                  <button type="button" className="task-item-body" onClick={() => onSelectJob(item)}>
                    <span className="task-item-top">
                      <strong>{item.type}</strong>
                      <code>{item.status}</code>
                    </span>
                    <span className="task-item-id">{item.id}</span>
                    <span className="task-item-time">{formatTime(item.updated_at)}</span>
                  </button>
                  <div className="task-item-actions">
                    <button type="button" className="ghost-mini" onClick={() => onSelectJob(item)}>
                      Open
                    </button>
                    <button
                      type="button"
                      className="danger-mini"
                      onClick={() => onDeleteJob(item.id)}
                      disabled={ACTIVE_JOB_STATUS.has(item.status)}
                      title={
                        ACTIVE_JOB_STATUS.has(item.status)
                          ? "Cancel running job before deleting"
                          : "Delete job"
                      }
                    >
                      Delete
                    </button>
                  </div>
                </div>
              ))
            ) : (
              <p className="hint">No tasks yet.</p>
            )}
          </div>
        </section>

        <section className="sidebar-block" aria-label="models">
          <h3>Models</h3>
          <div className="task-list">
            {recentModels.length > 0 ? (
              recentModels.map((item) => (
                <div key={item.id} className="task-item">
                  <div className="task-item-body task-item-body-static">
                    <span className="task-item-top">
                      <strong>model</strong>
                      <code>{item.id.slice(0, 8)}</code>
                    </span>
                    <span className="task-item-id">{item.path_ckpt}</span>
                    <span className="task-item-time">{formatTime(item.updated_at)}</span>
                  </div>
                  <div className="task-item-actions">
                    <button type="button" className="ghost-mini" onClick={() => openJobById(item.job_id)}>
                      Job
                    </button>
                    <a
                      className="ghost-mini anchor-mini"
                      href={buildApiUrl(`/api/results/models/${item.id}/download`)}
                    >
                      Download
                    </a>
                    <button type="button" className="danger-mini" onClick={() => onDeleteModel(item.id)}>
                      Delete
                    </button>
                  </div>
                </div>
              ))
            ) : (
              <p className="hint">No models yet.</p>
            )}
          </div>
        </section>

        <section className="sidebar-block" aria-label="results">
          <h3>Reports</h3>
          <div className="task-list">
            {recentResults.length > 0 ? (
              recentResults.map((item) => {
                const summaryAsset =
                  item.summary.assets.find((asset) => asset.name === "summary.json") ??
                  item.summary.assets[0] ??
                  null;
                return (
                  <div key={item.id} className="task-item">
                    <div className="task-item-body task-item-body-static">
                      <span className="task-item-top">
                        <strong>result</strong>
                        <code>{item.id.slice(0, 8)}</code>
                      </span>
                      <span className="task-item-id">{item.artifact_dir}</span>
                      <span className="task-item-time">{formatTime(item.updated_at)}</span>
                    </div>
                    <div className="task-item-actions">
                      <button
                        type="button"
                        className="ghost-mini"
                        onClick={() => openJobById(item.job_id)}
                      >
                        Job
                      </button>
                      {summaryAsset ? (
                        <a
                          className="ghost-mini anchor-mini"
                          href={buildApiUrl(summaryAsset.url)}
                        >
                          Download
                        </a>
                      ) : null}
                      <button type="button" className="danger-mini" onClick={() => onDeleteResult(item.id)}>
                        Delete
                      </button>
                    </div>
                  </div>
                );
              })
            ) : (
              <p className="hint">No reports yet.</p>
            )}
          </div>
        </section>
      </aside>
      <div className="workspace-main">
      <header className="hero">
        <h1>Squidiff LabFlow MVP</h1>
        <p>上传 → 校验 → Seurat 检查 → 500x500 预处理 → 训练任务 → 结果页</p>
      </header>

      <section className="panel">
        <h2>系统状态</h2>
        <p className="status">backend: {health}</p>
        {globalError ? <p className="error">{globalError}</p> : null}
      </section>

      {showTaskFeedback ? (
        <section className="panel task-feedback" aria-live="polite">
          <div className="task-progress" role="progressbar" aria-valuetext={taskLabel} />
          <p className="task-label">{taskLabel}</p>
        </section>
      ) : null}

      <section className="panel">
        <h2>1) 上传数据</h2>
        <div className="inline-guide-entry">
          <a
            className="guide-link-btn"
            href={buildApiUrl("/api/auth/user-guide")}
            target="_blank"
            rel="noreferrer"
          >
            User Guide
          </a>
        </div>
        <div className="form-grid">
          <label>
            数据名称 <ParamTooltip text={PARAM_TOOLTIPS["数据名称"]} />
            <input
              value={datasetName}
              onChange={(e) => setDatasetName(e.target.value)}
              placeholder="可选"
            />
          </label>
          <label>
            主数据文件（h5ad/rds/h5seurat） <ParamTooltip text={PARAM_TOOLTIPS["主数据文件"]} />
            <input
              type="file"
              accept=".h5ad,.rds,.h5seurat"
              onChange={(e) => setDatasetFile(e.target.files?.[0] ?? null)}
            />
          </label>
          <label>
            SMILES CSV（可选） <ParamTooltip text={PARAM_TOOLTIPS["SMILES CSV"]} />
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
            使用药物结构模式（需要 SMILES + dose 字段） <ParamTooltip text={PARAM_TOOLTIPS["使用药物结构模式"]} />
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
        <p className="hint">
          若 R 通过 Conda 安装（如 Windows），请选择 <strong>cmd_conda</strong>，并填写 conda.bat 的<strong>完整路径</strong>（如 F:\software\Miniconda3\condabin\conda.bat）和 R 环境名（如 r-4.3）。
        </p>
        <div className="form-grid compact">
          <label>
            R 执行方式 <ParamTooltip text={PARAM_TOOLTIPS["R 执行方式"]} />
            <select
              value={rExecMode}
              onChange={(e) => {
                setRExecMode(e.target.value as "direct" | "cmd_conda");
                setValidateError(null);
              }}
            >
              <option value="direct">direct (直接调用 Rscript)</option>
              <option value="cmd_conda">cmd_conda (cmd + conda activate)</option>
            </select>
          </label>
          <label>
            Rscript 命令 <ParamTooltip text={PARAM_TOOLTIPS["Rscript 命令"]} />
            <input
              value={rscriptBin}
              onChange={(e) => {
                setRscriptBin(e.target.value);
                setValidateError(null);
              }}
              placeholder="Rscript"
            />
          </label>
          <label>
            conda.bat 路径 <ParamTooltip text={PARAM_TOOLTIPS["conda.bat 路径"]} />
            {condaBatCandidates.length > 0 ? (
              <>
                <select
                  value={
                    condaBatCandidates.includes(rCondaBat) ? rCondaBat : "__custom__"
                  }
                  onChange={(e) => {
                    const v = e.target.value;
                    setValidateError(null);
                    if (v === "__custom__") {
                      setRCondaBat("");
                    } else {
                      setRCondaBat(v);
                      getCondaEnvs(v).then((res) =>
                        setCondaEnvsList(Array.from(new Set(res.conda_envs)))
                      );
                    }
                  }}
                >
                  {condaBatCandidates.map((path) => (
                    <option key={path} value={path}>
                      {path}
                    </option>
                  ))}
                  <option value="__custom__">其他（手动输入下方）</option>
                </select>
                {(!rCondaBat || !condaBatCandidates.includes(rCondaBat)) ? (
                  <input
                    value={rCondaBat}
                    onChange={(e) => {
                      const newBat = e.target.value;
                      setRCondaBat(newBat);
                      setValidateError(null);
                      // 如果输入了有效的路径，尝试获取环境列表
                      if (newBat.trim().length > 0 && newBat.endsWith("conda.bat")) {
                        getCondaEnvs(newBat.trim())
                          .then((res) =>
                            setCondaEnvsList(Array.from(new Set(res.conda_envs)))
                          )
                          .catch(() => {});
                      }
                    }}
                    placeholder="conda.bat 完整路径"
                  />
                ) : null}
              </>
            ) : (
              <input
                value={rCondaBat}
                onChange={(e) => {
                  setRCondaBat(e.target.value);
                  setValidateError(null);
                }}
                placeholder="conda.bat 或完整路径"
              />
            )}
          </label>
          <label>
            Conda R 环境名 <ParamTooltip text={PARAM_TOOLTIPS["Conda R 环境名"]} />
            {condaEnvsList.length > 0 ? (
              <select
                value={rCondaEnv}
                onChange={(e) => {
                  setRCondaEnv(e.target.value);
                  setValidateError(null);
                }}
              >
                <option value="">（请选择）</option>
                {[
                  ...(rCondaEnv && !condaEnvsList.includes(rCondaEnv)
                    ? [rCondaEnv]
                    : []),
                  ...condaEnvsList
                ].map((env) => (
                  <option key={env} value={env}>
                    {env}
                  </option>
                ))}
              </select>
            ) : (
              <input
                value={rCondaEnv}
                onChange={(e) => {
                  setRCondaEnv(e.target.value);
                  setValidateError(null);
                }}
                placeholder="例如: r-4.3"
              />
            )}
          </label>
        </div>
        <div className="step-actions">
          <button disabled={!canValidate || busyStep !== null} onClick={onValidate}>
            {busyStep === "validate" ? "校验中..." : "开始校验"}
          </button>
          {validateError ? (
            <div className="step-error" role="alert">
              <p className="step-error-msg">{validateError}</p>
              {validateRecommendation ? (
                <p className="step-error-recommend">{validateRecommendation}</p>
              ) : null}
            </div>
          ) : null}
        </div>
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
              <p className="hint">点击字段名可查看该列的分类取值，便于后续筛选细胞。</p>
              {seuratInspect.metadata_columns.length > 0 ? (
                <>
                  <div className="chip-list">
                    {seuratInspect.metadata_columns.map((column) => (
                      <button
                        key={column}
                        type="button"
                        className={`chip ${selectedMetadataColumn === column ? "chip-active" : ""}`}
                        onClick={() =>
                          setSelectedMetadataColumn((prev) =>
                            prev === column ? null : column
                          )
                        }
                        title={`查看 ${column} 的取值`}
                      >
                        {column}
                      </button>
                    ))}
                  </div>
                  {selectedMetadataColumn != null && (
                    <div className="metadata-values-panel">
                      <strong>{selectedMetadataColumn}</strong> 的取值（共{" "}
                      {(seuratInspect.metadata_column_values?.[selectedMetadataColumn]?.length ??
                        0)}{" "}
                      个分类）：
                      <div className="chip-list metadata-values-list">
                        {(seuratInspect.metadata_column_values?.[selectedMetadataColumn] ?? []).map(
                          (item) => (
                            <code key={item.value} className="value-chip" title={`${item.value}: ${item.count} 个细胞`}>
                              {item.value} <span className="value-count">({item.count})</span>
                            </code>
                          )
                        )}
                      </div>
                    </div>
                  )}
                </>
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
            分组字段（group_column） <ParamTooltip text={PARAM_TOOLTIPS["分组字段（group_column）"]} />
            <input
              value={groupColumn}
              onChange={(e) => setGroupColumn(e.target.value)}
              placeholder="例如: sample"
              list="metadata-columns"
            />
          </label>
          <label>
            类群字段（cluster_column） <ParamTooltip text={PARAM_TOOLTIPS["类群字段（cluster_column）"]} />
            <input
              value={clusterColumn}
              onChange={(e) => setClusterColumn(e.target.value)}
              placeholder="例如: celltype"
              list="metadata-columns"
            />
          </label>
          <label>
            随机种子（seed） <ParamTooltip text={PARAM_TOOLTIPS["随机种子（seed）"]} />
            <input
              type="number"
              value={prepareSeed}
              onChange={(e) => setPrepareSeed(Number(e.target.value))}
            />
          </label>
          <label>
            待筛选 clusters（逗号分隔） <ParamTooltip text={PARAM_TOOLTIPS["待筛选 clusters（逗号分隔）"]} />
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
            gene_size <ParamTooltip text={PARAM_TOOLTIPS.gene_size} />
            <input
              type="number"
              value={geneSize}
              onChange={(e) => setGeneSize(Number(e.target.value))}
            />
          </label>
          <label>
            output_dim <ParamTooltip text={PARAM_TOOLTIPS.output_dim} />
            <input
              type="number"
              value={outputDim}
              onChange={(e) => setOutputDim(Number(e.target.value))}
            />
          </label>
          <label>
            batch_size <ParamTooltip text={PARAM_TOOLTIPS.batch_size} />
            <input
              type="number"
              value={batchSize}
              onChange={(e) => setBatchSize(Number(e.target.value))}
            />
          </label>
          <label>
            lr <ParamTooltip text={PARAM_TOOLTIPS.lr} />
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
            {canCancelTrain ? (
              <div className="step-actions">
                <button type="button" className="danger-btn" onClick={onCancelTrain}>
                  {cancelBusy ? "停止中..." : "停止训练任务"}
                </button>
              </div>
            ) : null}
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
            {job.status === "running" ? (
              <section className="gpu-panel" aria-live="polite" aria-label="GPU stats">
                <div className="gpu-panel-head">
                  <h3>GPU Runtime</h3>
                  <span className="gpu-pulse" />
                  <span className="gpu-updated">
                    {gpuStats ? `updated ${formatTime(gpuStats.updated_at)}` : "polling..."}
                  </span>
                </div>
                {gpuStatsError ? <p className="error">{gpuStatsError}</p> : null}
                {gpuStats?.available === false ? (
                  <p className="hint">nvidia-smi unavailable: {gpuStats.reason ?? "unknown"}</p>
                ) : null}
                {gpuStats?.available && gpuStats.gpus.length > 0 ? (
                  <div className="gpu-grid">
                    {gpuStats.gpus.map((gpu) => {
                      const gpuUtil = clampPercent(gpu.utilization_gpu ?? 0);
                      const memPct =
                        gpu.memory_total_mb && gpu.memory_total_mb > 0
                          ? clampPercent((100 * (gpu.memory_used_mb ?? 0)) / gpu.memory_total_mb)
                          : 0;
                      return (
                        <article key={`${gpu.index ?? "gpu"}-${gpu.name}`} className="gpu-card">
                          <p className="gpu-title">
                            GPU {gpu.index ?? "-"}: {gpu.name}
                          </p>
                          <div className="gpu-row">
                            <span>Compute</span>
                            <strong>{Math.round(gpuUtil)}%</strong>
                          </div>
                          <div className="gpu-meter">
                            <span style={{ width: `${gpuUtil}%` }} />
                          </div>
                          <div className="gpu-row">
                            <span>VRAM</span>
                            <strong>
                              {gpu.memory_used_mb ?? 0}/{gpu.memory_total_mb ?? 0} MB
                            </strong>
                          </div>
                          <div className="gpu-meter gpu-meter-mem">
                            <span style={{ width: `${memPct}%` }} />
                          </div>
                          <div className="gpu-meta">
                            <span>mem util: {Math.round(gpu.utilization_memory ?? 0)}%</span>
                            <span>temp: {gpu.temperature_c ?? "-"}C</span>
                            <span>
                              power: {gpu.power_draw_w ?? "-"} / {gpu.power_limit_w ?? "-"} W
                            </span>
                          </div>
                        </article>
                      );
                    })}
                  </div>
                ) : (
                  <p className="hint">Waiting for GPU metrics...</p>
                )}
              </section>
            ) : null}
            {job.status === "canceled" ? (
              <p className="error">任务已由用户停止。</p>
            ) : null}
            <div className="live-log-box" aria-label="运行日志">
              <h3>运行日志</h3>
              <pre className="live-log-pre">
                {jobLog || liveJobLog || "（等待日志…）"}
              </pre>
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
      </div>
      </div>
      )}
    </main>
  );
}
