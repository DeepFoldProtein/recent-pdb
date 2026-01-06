#!/usr/bin/env python3
"""
In-House MSA Server
===================
A FastAPI-based MSA server compatible with ColabFold API format.
Uses MMseqs2 with optional GPU acceleration for high-performance MSA generation.

API Endpoints:
  POST /ticket/msa      - Submit MSA search job
  GET  /ticket/msa/{id} - Check job status / get results
  POST /result/msa      - Direct MSA search (synchronous)

Based on ColabFold MSA server protocol.
"""

import asyncio
import logging
import subprocess
import tempfile
import time
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

try:
    from fastapi import FastAPI, HTTPException, BackgroundTasks
    from fastapi.responses import PlainTextResponse
    from pydantic import BaseModel
except ImportError:
    print("FastAPI not installed. Run: pip install fastapi uvicorn")
    raise

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


# =============================================================================
# Configuration
# =============================================================================
@dataclass
class ServerConfig:
    """MSA Server configuration."""

    # MMseqs2 settings
    mmseqs_binary: str = "mmseqs"
    use_gpu: bool = False
    gpu_server_mode: bool = False

    # Database paths
    databases: dict = field(
        default_factory=lambda: {
            "uniref30": "/db/uniref30_2302",
            "colabfold_envdb": "/db/colabfold_envdb_202108",
        }
    )

    # Search parameters
    max_seqs: int = 10000
    sensitivity: float = 8.0

    # Server settings
    result_cache_dir: str = "/tmp/msa_cache"
    max_concurrent_jobs: int = 4
    job_timeout_seconds: int = 3600


config = ServerConfig()


# =============================================================================
# Job Management
# =============================================================================
@dataclass
class MSAJob:
    """MSA search job."""

    job_id: str
    sequence: str
    status: str = "pending"  # pending, running, complete, error
    created_at: float = field(default_factory=time.time)
    completed_at: Optional[float] = None
    result_path: Optional[str] = None
    error_message: Optional[str] = None


# In-memory job storage (use Redis for production)
jobs: dict[str, MSAJob] = {}
job_semaphore: Optional[asyncio.Semaphore] = None

# =============================================================================
# FastAPI App
# =============================================================================
app = FastAPI(
    title="PSP MSA Server",
    description="In-house MSA server compatible with ColabFold API",
    version="1.0.0",
)


class MSARequest(BaseModel):
    """MSA search request."""

    query: str  # FASTA format or raw sequence
    databases: list[str] = ["uniref30", "colabfold_envdb"]
    mode: str = "all"  # "all", "unpaired", "paired"


class MSATicketResponse(BaseModel):
    """Ticket response for async job."""

    id: str
    status: str


class MSAStatusResponse(BaseModel):
    """Job status response."""

    id: str
    status: str
    a3m: Optional[str] = None
    error: Optional[str] = None


# =============================================================================
# MMseqs2 Integration
# =============================================================================
def parse_fasta(fasta_str: str) -> list[tuple[str, str]]:
    """Parse FASTA string into list of (header, sequence) tuples."""
    sequences = []
    current_header = None
    current_seq = []

    for line in fasta_str.strip().split("\n"):
        line = line.strip()
        if line.startswith(">"):
            if current_header:
                sequences.append((current_header, "".join(current_seq)))
            current_header = line[1:]
            current_seq = []
        elif line:
            current_seq.append(line)

    if current_header:
        sequences.append((current_header, "".join(current_seq)))

    return sequences


def run_mmseqs_search(
    query_fasta: Path,
    output_a3m: Path,
    databases: list[str],
) -> bool:
    """Run MMseqs2 search and generate A3M output."""

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)

        # Run search for each database
        all_results = []

        for db_name in databases:
            if db_name not in config.databases:
                logger.warning(f"Unknown database: {db_name}")
                continue

            db_path = config.databases[db_name]
            result_prefix = tmp / f"result_{db_name}"

            # Build MMseqs2 command
            cmd = [
                config.mmseqs_binary,
                "easy-search",
                str(query_fasta),
                db_path,
                str(result_prefix) + ".m8",
                str(tmp / "tmp"),
                "-s",
                str(config.sensitivity),
                "--max-seqs",
                str(config.max_seqs),
                "--format-output",
                "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tseq",
            ]

            # Add GPU flags if enabled
            if config.use_gpu:
                cmd.extend(["--gpu", "1"])
                if config.gpu_server_mode:
                    cmd.extend(["--gpu-server", "1", "--db-load-mode", "2"])

            logger.info(f"Running: {' '.join(cmd)}")

            try:
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=config.job_timeout_seconds,
                )

                if result.returncode != 0:
                    logger.error(f"MMseqs2 error: {result.stderr}")
                    continue

                # Parse results
                m8_file = Path(str(result_prefix) + ".m8")
                if m8_file.exists():
                    with open(m8_file) as f:
                        for line in f:
                            parts = line.strip().split("\t")
                            if len(parts) >= 13:
                                all_results.append(
                                    {
                                        "target": parts[1],
                                        "pident": float(parts[2]),
                                        "evalue": float(parts[10]),
                                        "tseq": parts[12],
                                    }
                                )

            except subprocess.TimeoutExpired:
                logger.error(f"MMseqs2 timeout for {db_name}")
                continue

        # Convert to A3M format
        _write_a3m(query_fasta, all_results, output_a3m)

        return output_a3m.exists()


def _write_a3m(
    query_fasta: Path,
    results: list[dict],
    output_path: Path,
):
    """Write search results as A3M file."""
    # Read query sequence
    with open(query_fasta) as f:
        fasta_content = f.read()

    sequences = parse_fasta(fasta_content)
    if not sequences:
        return

    query_header, query_seq = sequences[0]

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        # Write query first
        f.write(f">{query_header}\n{query_seq}\n")

        # Write hits (sorted by identity, deduplicated)
        seen = {query_seq.upper().replace("-", "")}
        results_sorted = sorted(results, key=lambda x: -x["pident"])

        for i, hit in enumerate(results_sorted[: config.max_seqs]):
            seq = hit["tseq"]
            seq_key = seq.upper().replace("-", "")

            if seq_key not in seen:
                seen.add(seq_key)
                f.write(f">{hit['target']}\n{seq}\n")


async def process_msa_job(job: MSAJob, databases: list[str]):
    """Process MSA job asynchronously."""
    global job_semaphore

    if job_semaphore is None:
        job_semaphore = asyncio.Semaphore(config.max_concurrent_jobs)

    async with job_semaphore:
        job.status = "running"

        try:
            # Create temp files
            cache_dir = Path(config.result_cache_dir)
            cache_dir.mkdir(parents=True, exist_ok=True)

            query_fasta = cache_dir / f"{job.job_id}_query.fasta"
            output_a3m = cache_dir / f"{job.job_id}.a3m"

            # Write query
            sequences = parse_fasta(job.sequence)
            if not sequences:
                # Raw sequence without header
                sequences = [("query", job.sequence.replace("\n", "").replace(" ", ""))]

            with open(query_fasta, "w") as f:
                for header, seq in sequences:
                    f.write(f">{header}\n{seq}\n")

            # Run search in executor (blocking call)
            loop = asyncio.get_event_loop()
            success = await loop.run_in_executor(
                None,
                run_mmseqs_search,
                query_fasta,
                output_a3m,
                databases,
            )

            if success:
                job.status = "complete"
                job.result_path = str(output_a3m)
            else:
                job.status = "error"
                job.error_message = "MSA search failed"

            # Cleanup query file
            query_fasta.unlink(missing_ok=True)

        except Exception as e:
            logger.exception(f"Job {job.job_id} failed")
            job.status = "error"
            job.error_message = str(e)

        job.completed_at = time.time()


# =============================================================================
# API Endpoints
# =============================================================================
@app.post("/ticket/msa", response_model=MSATicketResponse)
async def create_msa_ticket(
    request: MSARequest,
    background_tasks: BackgroundTasks,
):
    """Submit MSA search job (async with ticket)."""
    job_id = str(uuid.uuid4())[:8]

    job = MSAJob(
        job_id=job_id,
        sequence=request.query,
    )
    jobs[job_id] = job

    # Start background task
    background_tasks.add_task(
        process_msa_job,
        job,
        request.databases,
    )

    return MSATicketResponse(id=job_id, status="pending")


@app.get("/ticket/msa/{job_id}", response_model=MSAStatusResponse)
async def get_msa_status(job_id: str):
    """Check MSA job status / get results."""
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")

    job = jobs[job_id]

    response = MSAStatusResponse(
        id=job.job_id,
        status=job.status,
        error=job.error_message,
    )

    # Include A3M result if complete
    if job.status == "complete" and job.result_path:
        result_path = Path(job.result_path)
        if result_path.exists():
            response.a3m = result_path.read_text()

    return response


@app.post("/result/msa", response_class=PlainTextResponse)
async def direct_msa_search(request: MSARequest):
    """Direct MSA search (synchronous, returns A3M)."""
    job_id = str(uuid.uuid4())[:8]

    job = MSAJob(
        job_id=job_id,
        sequence=request.query,
    )

    # Process synchronously
    await process_msa_job(job, request.databases)

    if job.status == "complete" and job.result_path:
        result_path = Path(job.result_path)
        if result_path.exists():
            return result_path.read_text()

    raise HTTPException(
        status_code=500, detail=job.error_message or "MSA search failed"
    )


@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {
        "status": "healthy",
        "gpu_enabled": config.use_gpu,
        "databases": list(config.databases.keys()),
        "active_jobs": sum(1 for j in jobs.values() if j.status == "running"),
    }


# =============================================================================
# CLI Entry Point
# =============================================================================
def main():
    import argparse

    parser = argparse.ArgumentParser(description="PSP MSA Server")
    parser.add_argument("--host", default="0.0.0.0", help="Bind host")
    parser.add_argument("--port", type=int, default=8000, help="Bind port")
    parser.add_argument("--mmseqs", default="mmseqs", help="MMseqs2 binary path")
    parser.add_argument("--gpu", action="store_true", help="Enable GPU acceleration")
    parser.add_argument("--gpu-server", action="store_true", help="Use GPU server mode")
    parser.add_argument("--db-uniref30", help="Path to UniRef30 database")
    parser.add_argument("--db-envdb", help="Path to ColabFold envdb")
    parser.add_argument(
        "--cache-dir", default="/tmp/msa_cache", help="Result cache directory"
    )
    parser.add_argument("--max-jobs", type=int, default=4, help="Max concurrent jobs")

    args = parser.parse_args()

    # Update config
    config.mmseqs_binary = args.mmseqs
    config.use_gpu = args.gpu
    config.gpu_server_mode = args.gpu_server
    config.result_cache_dir = args.cache_dir
    config.max_concurrent_jobs = args.max_jobs

    if args.db_uniref30:
        config.databases["uniref30"] = args.db_uniref30
    if args.db_envdb:
        config.databases["colabfold_envdb"] = args.db_envdb

    # Run server
    import uvicorn

    uvicorn.run(app, host=args.host, port=args.port)


if __name__ == "__main__":
    main()
