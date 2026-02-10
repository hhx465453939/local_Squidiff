from __future__ import annotations

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .api.datasets import router as datasets_router
from .api.jobs import router as jobs_router
from .api.results import router as results_router
from .api.seurat import router as seurat_router
from .runtime import job_queue

app = FastAPI(title="Squidiff LabFlow MVP", version="0.1.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(datasets_router)
app.include_router(jobs_router)
app.include_router(results_router)
app.include_router(seurat_router)


@app.on_event("startup")
def on_startup() -> None:
    job_queue.start()


@app.on_event("shutdown")
def on_shutdown() -> None:
    job_queue.stop()


@app.get("/api/health")
async def health() -> dict[str, str]:
    return {"status": "ok"}
