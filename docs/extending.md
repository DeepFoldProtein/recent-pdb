# Extension Guide

## Adding a New Model

### Step 1: Create Input Generator

Create `scripts/input_gen/mymodel.py`:

```python
from pathlib import Path
from .base import ModelInputGenerator, TargetFeatures, logger

class MyModelInputGenerator(ModelInputGenerator):
    name = "MyModel"
    input_format = "json"  # or pkl, npz
    
    def generate(self, target_features: TargetFeatures, output_dir: Path) -> Path:
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{target_features.target_id}.json"
        
        # Build input
        input_data = {
            "target_id": target_features.target_id,
            "chains": []
        }
        
        for chain in target_features.chains:
            chain_data = {"chain_id": chain.chain_id, "sequence": chain.sequence}
            if chain.msa_path:
                msa_seqs = self.load_msa(chain.msa_path)
                chain_data["msa"] = [seq for _, seq in msa_seqs]
            input_data["chains"].append(chain_data)
        
        import json
        with open(output_path, "w") as f:
            json.dump(input_data, f, indent=2)
        
        return output_path
```

### Step 2: Update Configuration

Add to `config.yaml`:

```yaml
enabled_models:
  - "MyModel"

models:
  MyModel:
    container: "${containers}/mymodel.sif"
    input_format: "json"
    gpu_required: true
```

### Step 3: Add Snakefile Rules

Add inference rule to `Snakefile`.

---

## Adding a New Scorer

### Step 1: Create Scorer Plugin

Create `scripts/eval/my_scorer.py`:

```python
from pathlib import Path
from .base import ScorerPlugin, EvalResult, logger

class MyScorer(ScorerPlugin):
    name = "my_score"
    description = "Custom quality metric"
    higher_is_better = True
    score_range = (0.0, 1.0)
    
    def score(self, prediction_path: Path, ground_truth_path: Path,
              target_id: str, model_name: str) -> EvalResult:
        
        if not self.validate_structures(prediction_path, ground_truth_path):
            return EvalResult(
                scorer_name=self.name,
                target_id=target_id,
                model_name=model_name,
                score=0.0,
                metadata={"error": "Validation failed"}
            )
        
        # Your scoring logic here
        score = compute_my_metric(prediction_path, ground_truth_path)
        
        return EvalResult(
            scorer_name=self.name,
            target_id=target_id,
            model_name=model_name,
            score=score
        )
```

### Step 2: Update Configuration

```yaml
enabled_scorers:
  - "my_score"
```

---

## Using OpenStructure via Apptainer

The lDDT scorer uses OpenStructure container:

```bash
# Build container
apptainer build ost.sif \
    docker://registry.scicore.unibas.ch/schwede/openstructure:latest

# Run lDDT calculation
apptainer exec ost.sif ost scripts/eval/ost_scorer.py \
    --prediction pred.cif \
    --ground-truth ref.cif.gz \
    --output result.json
```

---

## Plugin Discovery

Plugins are discovered by naming convention:

- Input generators: `scripts/input_gen/{model_name.lower()}.py`
- Scorers: `scripts/eval/{scorer_name}_scorer.py`

---

## Testing Plugins

```bash
# Unit test
uv run pytest tests/test_my_scorer.py -v

# Integration test
uv run snakemake --dry-run \
    outputs/evaluation/AlphaFold3/8WGV/my_score.json
```
