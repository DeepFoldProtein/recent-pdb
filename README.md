# PSP Benchmark Pipeline

High-Performance Protein Structure Prediction Benchmark

## Quick Start

```bash
uv sync
uv run python scripts/setup_test.py
uv run snakemake -n
```

## Models

- AlphaFold3 (Apptainer)
- Protenix (Apptainer)
- ColabFold (Apptainer)

## Scorers

- lDDT (OpenStructure via Apptainer)
- DockQ
- TM-score
