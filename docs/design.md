# System Design Document

## Overview

PSP Benchmark Pipeline is a high-performance system for evaluating protein structure prediction models.

## Design Goals

1. **Reproducibility**: Identical inputs produce identical results
2. **Efficiency**: Sequence-based MSA caching prevents redundant computation
3. **Extensibility**: Plugin architecture for models and scorers
4. **Idempotency**: Resume from interruption without re-running completed tasks

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    Snakemake Orchestrator                   │
├─────────────────────────────────────────────────────────────┤
│  Stage 1   │ Stage 2  │ Stage 3 │ Stage 4 │   Stage 5       │
│  PDB Sync  │ Registry │   MSA   │  Input  │  Inference      │
│            │          │ Search  │   Gen   │  (Apptainer)    │
├────────────┴──────────┴─────────┴─────────┴─────────────────┤
│                     Shared MSA Cache                        │
│              (Sequence Hash → MSA Results)                  │
├─────────────────────────────────────────────────────────────┤
│  Stage 6: Evaluation          │  Stage 7: Aggregation       │
│  ┌─────┐ ┌─────┐ ┌─────┐     │  ┌─────────────────────┐     │
│  │lDDT │ │DockQ│ │ TM  │     │  │ Leaderboard (CSV)   │     │
│  │(OST)│ │     │ │score│     │  │              (HTML) │     │
│  └─────┘ └─────┘ └─────┘     │  │              (JSON) │     │
└─────────────────────────────────────────────────────────────┘
         ↑
    ┌────┴───────┐
    │ MSA Server │  (In-house, GPU-accelerated)
    └────────────┘
```

## Key Components

### 1. Sequence-Based Hashing

All protein sequences are identified by SHA256 hash:

```python
seq_hash = sha256(sequence.upper()).hexdigest()[:16]
```

**Benefits:**

- Search MSA only once per unique sequence
- Automatic handling of homo-dimers
- O(1) cache lookup

### 2. Model Integration (Apptainer)

Models run in isolated containers:

| Model | Container | Input | Output |
|-------|-----------|-------|--------|
| AlphaFold3 | `alphafold3.sif` | JSON | CIF |
| Protenix | `protenix.sif` | JSON | CIF |
| Protenix | `protenix.sif` | JSON | CIF |
| ColabFold | `colabfold.sif` | A3M/FASTA | PDB |

### 3. Featurizers (New)

Flexible MSA generation strategy:

- Combine multiple aligners (MMseqs2, HHblits)
- Combine multiple databases (UniRef30, PDB100)
- Configurable merge strategies (Top-N, Uniform, etc.)

### 4. In-House MSA Server

FastAPI-based server compatible with ColabFold API:

- **Endpoints**: `/ticket/msa`, `/result/msa`, `/health`
- **GPU Support**: MMseqs2 GPU server mode
- **Databases**: UniRef30, ColabFold EnvDB

### 5. Evaluation Scorers

| Scorer | Metric | Tool |
|--------|--------|------|
| OST | lDDT | OpenStructure (Apptainer) |
| DockQ | Interface Quality | DockQ.py |
| TM-score | Global Similarity | TMalign |

### 6. Snakemake DAG

```
sync_pdb → build_registry → run_msa_* → merge_msa → gen_*_input
                                                         ↓
                      summary_report ← run_scorer ← run_*
```

## Data Flow

```
PDB Mirror (.cif.gz)
    ↓
[sync_pdb.py] → SQLite Registry
    ↓
[target_registry.py] → targets.txt, seq_hashes.txt, sequences/
    ↓
[msa_search.py] → cache/{hash}/msa/merged.a3m
    ↓
[input_gen/*.py] → outputs/input/{model}/{target}/
    ↓
[Apptainer Container] → outputs/predictions/{model}/{target}.cif
    ↓
[eval/*.py] → outputs/evaluation/{model}/{target}/{scorer}.json
    ↓
[aggregate.py] → outputs/leaderboard.{csv,html,json}
```

## Container Registry

| Image | Source |
|-------|--------|
| AlphaFold3 | `docker://google-deepmind/alphafold3` |
| Protenix | `docker://bytedance/protenix` |
| ColabFold | `docker://ghcr.io/sokrypton/colabfold` |
| OpenStructure | `registry.scicore.unibas.ch/schwede/openstructure` |
| MMseqs2 | `docker://soedinglab/mmseqs2` |
| HH-suite | `docker://soedinglab/hh-suite` |

## Error Handling

| Error | Strategy |
|-------|----------|
| MSA Timeout | Log, skip, use empty MSA |
| Inference OOM | Retry with smaller batch |
| Scorer Failure | Record 0, save error message |
| SLURM Timeout | `--rerun-incomplete` |
