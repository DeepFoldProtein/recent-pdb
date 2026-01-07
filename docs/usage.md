# Usage Guide

## Prerequisites

- Python ≥ 3.11
- [uv](https://github.com/astral-sh/uv)
- Docker (for inference containers)
- Apptainer/Singularity
- SLURM (cluster execution)

## Installation

```bash
git clone <repo-url>
cd psp-benchmark
uv sync
```

## Environment Setup

Copy and configure the environment file:

```bash
cp setup_env.sh.template setup_env.sh
# Edit setup_env.sh with your paths:
#   PDB_MMCIF_DIR, AF3_DB_DIR, AF3_MODEL_DIR, etc.

source setup_env.sh
```

## Quick Test

```bash
# Setup test data
uv run python scripts/setup_test.py

# Dry-run
uv run snakemake -n
```

---

## Running the Full Pipeline

### Single Target (End-to-End)

Run the complete pipeline for a single target (e.g., `4HHB`):

```bash
# Set environment
source setup_env.sh

# Run MSA → Input Gen → Inference → Evaluation for AlphaFold3
uv run snakemake --cores 8 --use-singularity \
    outputs/evaluation/AlphaFold3/4HHB/ost.json
```

This will automatically:

1. Generate MSAs (MMseqs2 uniref30 + pdb100)
2. Create AlphaFold3 input JSON
3. Run AF3 inference
4. Evaluate with OST (lDDT, QS-score, DockQ, etc.)

### Individual Steps

```bash
# MSA only
uv run snakemake --cores 8 cache/sequences/mmseqs2_18/uniref30_2023_02/.../{seq_hash}.a3m

# Input generation only
uv run snakemake --cores 8 outputs/input/AlphaFold3/{target_id}/{target_id}.json

# Inference only (requires input)
uv run snakemake --cores 8 --use-singularity outputs/predictions/AlphaFold3/{target_id}.cif

# Evaluation only (requires prediction)
uv run snakemake --cores 8 --use-singularity outputs/evaluation/AlphaFold3/{target_id}/ost.json
```

### All Targets

```bash
# Run evaluation for all enabled models and targets
uv run snakemake run_evaluation --cores 8 --use-singularity
```

## Container Building

### MSA & Evaluation Containers (Direct Pull)

These containers can be built directly without Docker:

```bash
./scripts/build_containers.sh mmseqs2 hhsuite hmmer openstructure
```

### Inference Containers (Require Docker)

AlphaFold3, Protenix, ColabFold require Docker build then Apptainer conversion:

```bash
# Requires Docker with GPU support
./scripts/build_docker_containers.sh alphafold3
./scripts/build_docker_containers.sh colabfold
./scripts/build_docker_containers.sh protenix

# Note: Protenix requires a one-time resource setup (CCD, etc.)
./scripts/setup_resources.sh
```

**Container Status:**

| Container | Method | Source |
| :--- | :--- | :--- |
| mmseqs2 | Direct pull | ghcr.io/soedinglab |
| hhsuite | Direct pull | soedinglab/hh-suite |
| hmmer | Direct pull | biocontainers/hmmer |
| openstructure | Direct pull | registry.scicore.unibas.ch |
| alphafold3 | Docker build | github.com/google-deepmind/alphafold3 |
| colabfold | Docker build | github.com/sokrypton/ColabFold |
| protenix | Docker build | github.com/bytedance/Protenix |

### Verification

After building, verify that all containers are functional:

```bash
./scripts/verify_containers.sh
```

---

## Multi-Cluster Execution

### CPU Cluster (MSA, Evaluation)

```bash
uv run snakemake run_all_msa --profile slurm/ --jobs 100
uv run snakemake run_evaluation --profile slurm/ --jobs 50
```

### GPU Cluster (Inference)

```bash
uv run snakemake run_alphafold3 run_protenix run_colabfold \
    --profile slurm-gpu/ --jobs 10
```

Create `slurm-gpu/config.yaml`:

```yaml
executor: slurm
default-resources:
  slurm_partition: gpu
  slurm_extra: "'--gres=gpu:1'"
```

---

## Pipeline Stages

| Stage | Rule | Cluster |
| :--- | :--- | :--- |
| 1-2 | `sync_pdb`, `build_registry` | Any |
| 3 | `run_all_msa` | CPU (high mem) |
| 4 | `gen_*_input` | CPU |
| 5 | `run_alphafold3`, etc. | GPU |
| 6-7 | `run_evaluation`, `summary_report` | CPU |

---

## Custom Target Lists

Instead of automatic PDB filtering, you can provide a custom list of PDB IDs.

### Option 1: Manual List

Create a text file with one PDB ID per line:

```text
# my_targets.txt
# Comments start with #
8abc
9xyz
7def
```

Then set in `config.yaml`:

```yaml
filters:
  custom_targets_file: "./my_targets.txt"
```

### Option 2: Generate with Script

Use `generate_pdb_ids.py` to create lists based on criteria:

```bash
# Sample 100 PDBs from 2024 with resolution ≤ 2.0
uv run python scripts/generate_pdb_ids.py \
    --registry-db data/registry.sqlite \
    --date-range 2024-01-01 2024-12-31 \
    --resolution-max 2.0 \
    --sample 100 \
    --output my_targets.txt

# Exclude specific PDBs
uv run python scripts/generate_pdb_ids.py \
    --registry-db data/registry.sqlite \
    --exclude obsolete.txt \
    --sample 50 \
    --output filtered.txt
```

**Options:**

| Flag | Description |
| :--- | :--- |
| `--resolution-max R` | Maximum resolution in Å |
| `--date-range START END` | Release date range (YYYY-MM-DD) |
| `--sample N` | Random sample N entries |
| `--exclude FILE` | Exclude PDB IDs from file |
| `--seed N` | Random seed for reproducibility |

---

## Configuration

Edit `config.yaml`:

```yaml
paths:
  msa_cache: "./cache/sequences"

msa:
  tools:
    mmseqs2:
      databases:
        - name: "uniref50"
          path: "/store/database/mmseqs/uniref50"

enabled_models:
  - "AlphaFold3"
  - "Protenix"
  - "ColabFold"

models:
  AlphaFold3:
    use_templates: true # Set to false to disable template generation

# -----------------------------------------------------------------------------
# Featurizers (New)
# -----------------------------------------------------------------------------
# Define how to search and merge MSAs
featurizers:
  default_setup:
    sources:
      - aligner: mmseqs2
        database: uniref30_2023_02
    merge:
      strategy: "top_n"

enabled_featurizers:
  - "default_setup"
```

---

## PDB Mirror Management

The pipeline includes a built-in mechanism to synchronize and index the entire PDB archive.

### Full Synchronization

To download (or update) the entire PDB mmCIF archive and rebuild the metadata registry:

```bash
# Run full rsync and registry update
uv run python scripts/sync_pdb.py \
    --pdb-master ${PDB_MMCIF_DIR} \
    --registry-db data/registry.sqlite
```

> [!NOTE]
> The default `Snakefile` rule `sync_pdb` uses the `--skip-sync` flag to save time during routine execution. Use the command above for periodic updates.

### Force Rescan

If you suspect the registry is out of sync with the files on disk, or after modifying the schema (e.g., adding SEQRES):

```bash
uv run python scripts/sync_pdb.py \
    --registry-db data/registry.sqlite \
    --skip-sync \
    --force-rescan
```

## Data Management

The pipeline distinguishes between different types of data to prevent accidental deletion of expensive resources.

### Directory Structure

- `outputs/`: Volatile results (predictions, reports). Safe to delete.
- `cache/sequences/`: Expensive-to-regenerate MSA files.
- `data/registry.sqlite`: The metadata registry.
- `containers/*.sif`: Apptainer ready-to-use images.
- `resources/`: Persistent resources (model weights, CCD files) managed by `scripts/setup_resources.sh`.

### Cleaning Rules

Use these Snakemake rules to clean specific parts of the pipeline:

```bash
# Remove only predictions and evaluation reports (safe)
uv run snakemake clean_outputs

# Remove MSA cache (will require rerunning all MSA searches)
uv run snakemake clean_msa

# Remove registry and target lists (will require rerunning sync_pdb)
uv run snakemake clean_registry

# [DANGEROUS] Remove all generated data (outputs and caches)
uv run snakemake clean_all
```

> [!IMPORTANT]
> None of the `clean` rules will delete your `containers/` or the resources in `resources/` (model weights/CCD). These must be managed manually or via the management script:
>
> ```bash
> ./scripts/manage_resources.sh --status
> ```
