# =============================================================================
# PSP Benchmark Pipeline - Main Snakemake Workflow
# =============================================================================
# High-performance protein structure prediction benchmark pipeline.
# Orchestrates 7 stages from PDB sync to final leaderboard generation.
#
# Usage:
#   snakemake --cores 8                    # Local execution
#   snakemake --profile slurm/             # SLURM cluster execution
#   snakemake -n --printshellcmds          # Dry-run with commands
#
# Apptainer Usage:
#   snakemake --use-singularity --cores 8  # Use Apptainer containers
# =============================================================================

import os
import json
from pathlib import Path

def deep_merge(base, override):
    """Recursively merge override dict into base dict."""
    for key, value in override.items():
        if key in base and isinstance(base[key], dict) and isinstance(value, dict):
            deep_merge(base[key], value)
        else:
            base[key] = value
    return base

# Load configuration
configfile: "config.yaml"

# Merge local overrides if exists
if os.path.exists("config.local.yaml"):
    import yaml
    with open("config.local.yaml") as f:
        local_config = yaml.safe_load(f) or {}
    deep_merge(config, local_config)

# Resolve paths (relative for portability)
CACHE_DIR = Path(config["paths"]["msa_cache"])
OUTPUT_DIR = Path(config["paths"]["output_dir"])
REGISTRY_DB = Path(config["paths"]["registry_db"])
TARGET_LISTS = Path(config["paths"]["target_lists"])
GROUND_TRUTHS = Path(config["paths"]["ground_truths"])
CONTAINERS = Path(config["paths"].get("containers", "./containers"))

# Enabled models, scorers, and featurizers
ENABLED_MODELS = config.get("enabled_models", ["AlphaFold3", "Protenix", "ColabFold"])
ENABLED_SCORERS = config.get("enabled_scorers", ["ost", "dockq", "tm_score"])
ENABLED_FEATURIZERS = config.get("enabled_featurizers", [])

# =============================================================================
# Dynamic Target Discovery
# =============================================================================
def get_target_ids():
    """Load target IDs from registry."""
    targets_file = TARGET_LISTS / "targets.txt"
    if targets_file.exists():
        return [line.strip() for line in targets_file.open() if line.strip()]
    return []

def get_seq_hashes():
    """Load unique sequence hashes from registry."""
    hashes_file = TARGET_LISTS / "seq_hashes.txt"
    if hashes_file.exists():
        return [line.strip() for line in hashes_file.open() if line.strip()]
    return []

def get_target_seq_hashes(target_id):
    """Get sequence hashes for a target (from mapping file)."""
    mapping_file = TARGET_LISTS / "seq_to_targets.json"
    if mapping_file.exists():
        mapping = json.loads(mapping_file.read_text())
        return [h for h, targets in mapping.items() if target_id in targets]
    return []

# =============================================================================
# Container Definitions (Apptainer/Singularity)
# =============================================================================
container: None  # Default: no container

# =============================================================================
# Stage 1: PDB Sync & Audit
# =============================================================================
rule sync_pdb:
    """Synchronize PDB mirror and update metadata database."""
    output:
        touch(OUTPUT_DIR / ".sync_done")
    resources:
        cpus=1,
        mem_mb=4000,
        runtime=120
    shell:
        """
        python scripts/sync_pdb.py \
            --registry-db {REGISTRY_DB} \
            --config config.yaml \
            --skip-sync
        """

# =============================================================================
# Stage 2: Target Registry
# =============================================================================
rule build_registry:
    """Filter targets and generate sequence hash mappings."""
    input:
        rules.sync_pdb.output
    output:
        targets=TARGET_LISTS / "targets.txt",
        hashes=TARGET_LISTS / "seq_hashes.txt",
        mapping=TARGET_LISTS / "seq_to_targets.json"
    resources:
        cpus=2,
        mem_mb=8000,
        runtime=30
    shell:
        """
        python scripts/target_registry.py \
            --registry-db {REGISTRY_DB} \
            --config config.yaml \
            --output-dir {TARGET_LISTS}
        """

# =============================================================================
# Stage 3: MSA Search (MMseqs2 + HHblits)
# =============================================================================

def get_db_map(tool_name):
    """Map directory name -> db config."""
    mapping = {}
    for db in config["msa"]["tools"].get(tool_name, {}).get("databases", []):
        v = db.get("version", "")
        key = f"{db['name']}_{v}" if v else db['name']
        mapping[key] = db
    return mapping

def get_db_path(wildcards, tool_name):
    return get_db_map(tool_name)[wildcards.db_dir]["path"]

def get_db_args(wildcards, tool_name):
    db = get_db_map(tool_name)[wildcards.db_dir]
    return db["name"], db.get("version", "")

rule run_msa_mmseqs2:
    """Run MMseqs2 search for a single sequence hash and database."""
    input:
        query=TARGET_LISTS / "sequences" / "{seq_hash}.fasta"
    output:
        a3m=CACHE_DIR / "mmseqs2_{version}" / "{db_dir}" / "{phash}" / "{seq_hash}.a3m"
    params:
        cache_dir=CACHE_DIR,
        container=CONTAINERS / "mmseqs2.sif" if (CONTAINERS / "mmseqs2.sif").exists() else None,
        db_path=lambda wc: get_db_path(wc, "mmseqs2"),
        msa_params=lambda wc: config["msa"]["tools"]["mmseqs2"]["params"]
    resources:
        cpus=16,
        mem_mb=32000,
        runtime=240
    run:
        db_name, db_version = get_db_args(wildcards, "mmseqs2")
        container_arg = f"--container {params.container}" if params.container else ""
        shell(f"""
        python -m scripts.msa_search \
            --seq-hash {wildcards.seq_hash} \
            --cache-dir {params.cache_dir} \
            --config config.yaml \
            --tool mmseqs2 \
            --db-name {db_name} \
            --db-version "{db_version}" \
            {container_arg} \
            --threads {resources.cpus} \
            --output {output.a3m}
        """)

rule run_msa_hhblits:
    """Run HHblits search for a single sequence hash and database."""
    input:
        query=TARGET_LISTS / "sequences" / "{seq_hash}.fasta"
    output:
        a3m=CACHE_DIR / "hhblits_{version}" / "{db_dir}" / "{phash}" / "{seq_hash}.a3m"
    params:
        cache_dir=CACHE_DIR,
        container=CONTAINERS / "hhsuite.sif" if (CONTAINERS / "hhsuite.sif").exists() else None,
        db_path=lambda wc: get_db_path(wc, "hhblits"),
        msa_params=lambda wc: config["msa"]["tools"]["hhblits"]["params"]
    resources:
        cpus=16,
        mem_mb=32000,
        runtime=240
    run:
        db_name, db_version = get_db_args(wildcards, "hhblits")
        shell(f"""
        mkdir -p $(dirname {output.a3m})
        python -m scripts.msa_search \
            --seq-hash {wildcards.seq_hash} \
            --cache-dir {params.cache_dir} \
            --config config.yaml \
            --tool hhblits \
            --db-name {db_name} \
            --db-version "{db_version}" \
            --threads {resources.cpus} \
            --output {output.a3m}
        """)

def get_params_hash(tool_name):
    """Compute short hash of tool parameters."""
    params = config["msa"]["tools"].get(tool_name, {}).get("params", {})
    if params is None:
        params = {}
    s = json.dumps(params, sort_keys=True)
    import hashlib
    return hashlib.md5(s.encode()).hexdigest()[:8]

def get_tool_version(tool_name):
    return config["msa"]["tools"].get(tool_name, {}).get("version", "unknown")

def get_msa_files_for_hash(seq_hash):
    """Get all expected MSA files for a sequence hash."""
    files = []
    
    # helper to build path
    def _add_files(tool_name):
        version = get_tool_version(tool_name)
        phash = get_params_hash(tool_name)
        tool_dir = f"{tool_name}_{version}"
        for key in get_db_map(tool_name):
             files.append(CACHE_DIR / tool_dir / key / phash / f"{seq_hash}.a3m")

    for tool in config["msa"]["tools"]:
        if config["msa"]["tools"][tool].get("enabled", True):
            _add_files(tool)
    
    return files

def get_target_msa_files(target_id):
    """Get all MSA files for all sequences in a target."""
    files = []
    for h in get_target_seq_hashes(target_id):
        files.extend(get_msa_files_for_hash(h))
    return files

def get_featurizer_msa_files(featurizer_name, seq_hash):
    """Get MSA files for a specific featurizer's sources."""
    featurizer_cfg = config.get("featurizers", {}).get(featurizer_name, {})
    sources = featurizer_cfg.get("sources", [])
    files = []
    for src in sources:
        aligner = src["aligner"]
        db_name = src["database"]
        version = get_tool_version(aligner)
        phash = get_params_hash(aligner)
        # Build db_dir key (name_version format)
        db_map = get_db_map(aligner)
        for db_key, db_cfg in db_map.items():
            if f"{db_cfg['name']}_{db_cfg.get('version', '')}" == db_name or db_key == db_name:
                files.append(CACHE_DIR / f"{aligner}_{version}" / db_key / phash / f"{seq_hash}.a3m")
                break
    return files

def get_target_featurizer_msa_files(featurizer_name, target_id):
    """Get all MSA files for a featurizer across all sequences in a target."""
    files = []
    for h in get_target_seq_hashes(target_id):
        files.extend(get_featurizer_msa_files(featurizer_name, h))
    return files

def get_all_msa_files():
    """Get all MSA files for all sequence hashes."""
    files = []
    for seq_hash in get_seq_hashes():
        files.extend(get_msa_files_for_hash(seq_hash))
    return files

rule run_all_msa:
    """Aggregate rule for all MSA search jobs."""
    input:
        get_all_msa_files()

# =============================================================================
# Stage 4: Model Input Generation
# =============================================================================
rule gen_alphafold3_input:
    """Generate AlphaFold3 JSON input."""
    input:
        msa_done=lambda wc: get_target_featurizer_msa_files(wc.featurizer, wc.target_id)
    output:
        json=OUTPUT_DIR / "input" / "{featurizer}" / "AlphaFold3" / "{target_id}" / "{target_id}.json"
    shell:
        """
        python -m scripts.input_gen.alphafold3 \
            --target-id {wildcards.target_id} \
            --featurizer {wildcards.featurizer} \
            --cache-dir {CACHE_DIR} \
            --output-dir {OUTPUT_DIR}/input/{wildcards.featurizer}/AlphaFold3 \
            --config config.yaml
        """

rule gen_protenix_input:
    """Generate Protenix JSON input."""
    input:
        msa_done=lambda wc: get_target_featurizer_msa_files(wc.featurizer, wc.target_id)
    output:
        json=OUTPUT_DIR / "input" / "{featurizer}" / "Protenix" / "{target_id}" / "{target_id}.json"
    shell:
        """
        python -m scripts.input_gen.protenix \
            --target-id {wildcards.target_id} \
            --featurizer {wildcards.featurizer} \
            --cache-dir {CACHE_DIR} \
            --output-dir {OUTPUT_DIR}/input/{wildcards.featurizer}/Protenix \
            --config config.yaml
        """

rule gen_colabfold_input:
    """Generate ColabFold FASTA/A3M input."""
    input:
        msa_done=lambda wc: get_target_featurizer_msa_files(wc.featurizer, wc.target_id)
    output:
        fasta=OUTPUT_DIR / "input" / "{featurizer}" / "ColabFold" / "{target_id}" / "{target_id}.a3m"
    shell:
        """
        python -m scripts.input_gen.colabfold \
            --target-id {wildcards.target_id} \
            --featurizer {wildcards.featurizer} \
            --cache-dir {CACHE_DIR} \
            --output-dir {OUTPUT_DIR}/input/{wildcards.featurizer}/ColabFold \
            --config config.yaml
        """

# =============================================================================
# Stage 5: Inference (Apptainer Containers)
# =============================================================================
rule run_alphafold3:
    """Run AlphaFold3 inference."""
    input:
        json=OUTPUT_DIR / "input" / "{featurizer}" / "AlphaFold3" / "{target_id}" / "{target_id}.json"
    output:
        cif=OUTPUT_DIR / "pred" / "{featurizer}" / "AlphaFold3" / "{target_id}.cif"
    params:
        container=lambda wc: config["models"]["AlphaFold3"]["container"],
        model_dir=lambda wc: config["models"]["AlphaFold3"]["model_dir"],
        db_dir=lambda wc: config["models"]["AlphaFold3"]["databases_dir"]
    resources:
        cpus=4,
        mem_mb=80000,
        runtime=120,
        gpu=1
    shell:
        """
        mkdir -p $(dirname {output.cif})
        
        apptainer exec --nv \
            --bind $(dirname {input.json}):/input \
            --bind $(dirname {output.cif}):/output \
            --bind {params.model_dir}:/root/models \
            --bind {params.db_dir}:/root/public_databases \
            {params.container} \
            python /app/alphafold/run_alphafold.py \
                --json_path=/input/$(basename {input.json}) \
                --model_dir=/root/models \
                --output_dir=/output \
                --run_data_pipeline=false \
                --run_inference=true
        
        # Rename output to expected path (handles timestamped dirs)
        SOURCE_CIF=$(find $(dirname {output.cif}) -name "{wildcards.target_id}_model.cif" -printf "%T@ %p\n" | sort -n | tail -1 | cut -f2- -d" ")
        if [ ! -z "$SOURCE_CIF" ]; then
            mv "$SOURCE_CIF" {output.cif}
        fi
        """

rule run_protenix:
    """Run Protenix inference."""
    input:
        json=OUTPUT_DIR / "input" / "{featurizer}" / "Protenix" / "{target_id}" / "{target_id}.json"
    output:
        cif=OUTPUT_DIR / "pred" / "{featurizer}" / "Protenix" / "{target_id}.cif"
    params:
        container=lambda wc: config["models"]["Protenix"]["container"],
        model_name=lambda wc: config["models"]["Protenix"]["model_name"]
    resources:
        cpus=4,
        mem_mb=80000,
        runtime=120,
        gpu=1
    shell:
        """
        mkdir -p $(dirname {output.cif})
        
        apptainer exec --nv \
            --bind $(dirname {input.json}):/input \
            --bind $(dirname {output.cif}):/output \
            {params.container} \
            protenix predict \
                --input /input/$(basename {input.json}) \
                --out_dir /output \
                --model_name {params.model_name} \
                --seeds 101
        
        # Rename output
        mv $(dirname {output.cif})/{wildcards.target_id}/*.cif {output.cif} 2>/dev/null || \
        mv $(dirname {output.cif})/*.cif {output.cif} 2>/dev/null || true
        """

rule run_colabfold:
    """Run ColabFold inference."""
    input:
        a3m=OUTPUT_DIR / "input" / "{featurizer}" / "ColabFold" / "{target_id}" / "{target_id}.a3m"
    output:
        pdb=OUTPUT_DIR / "pred" / "{featurizer}" / "ColabFold" / "{target_id}.cif"
    params:
        container=lambda wc: config["models"]["ColabFold"]["container"],
        model_type=lambda wc: config["models"]["ColabFold"]["model_type"],
        num_recycle=lambda wc: config["models"]["ColabFold"]["num_recycle"]
    resources:
        cpus=4,
        mem_mb=32000,
        runtime=120,
        gpu=1
    shell:
        """
        mkdir -p $(dirname {output.pdb})
        
        apptainer exec --nv \
            --bind $(dirname {input.a3m}):/input \
            --bind $(dirname {output.pdb}):/output \
            {params.container} \
            colabfold_batch \
                /input/$(basename {input.a3m}) \
                /output \
                --model-type {params.model_type} \
                --num-recycle {params.num_recycle}
        
        # Convert PDB to CIF if needed, or rename
        mv $(dirname {output.pdb})/{wildcards.target_id}/*_relaxed_rank_001*.pdb {output.pdb} 2>/dev/null || \
        mv $(dirname {output.pdb})/*_relaxed_rank_001*.pdb {output.pdb} 2>/dev/null || \
        mv $(dirname {output.pdb})/*_unrelaxed_rank_001*.pdb {output.pdb} 2>/dev/null || true
        """

rule run_all_inference:
    """Aggregate rule for all inference jobs."""
    input:
        expand(OUTPUT_DIR / "pred" / "{featurizer}" / "{model}" / "{target_id}.cif",
               featurizer=ENABLED_FEATURIZERS,
               model=ENABLED_MODELS,
               target_id=get_target_ids())

# =============================================================================
# Stage 6: Evaluation (Scorer Plugins)
# =============================================================================

def get_ground_truth(wildcards):
    """Resolve ground truth path from registry mapping."""
    mapping_file = TARGET_LISTS / "ground_truths.json"
    if mapping_file.exists():
        mapping = json.loads(mapping_file.read_text())
        return mapping.get(wildcards.target_id.upper())
    return ""

rule run_scorer:
    """Run a single scorer on a prediction."""
    input:
        prediction=OUTPUT_DIR / "pred" / "{featurizer}" / "{model}" / "{target_id}.cif",
        ground_truth=get_ground_truth
    output:
        score=OUTPUT_DIR / "eval" / "{featurizer}" / "{model}" / "{target_id}" / "{scorer}.json"
    resources:
        cpus=2,
        mem_mb=8000,
        runtime=30
    shell:
        """
        uv run python -m scripts.eval.{wildcards.scorer}_scorer \
            --prediction {input.prediction} \
            --ground-truth {input.ground_truth} \
            --output {output.score}
        """

rule run_evaluation:
    """Aggregate rule for all evaluation jobs."""
    input:
        expand(OUTPUT_DIR / "eval" / "{featurizer}" / "{model}" / "{target_id}" / "{scorer}.json",
               featurizer=ENABLED_FEATURIZERS,
               model=ENABLED_MODELS,
               target_id=get_target_ids(),
               scorer=ENABLED_SCORERS)

# =============================================================================
# Stage 7: Aggregation & Leaderboard
# =============================================================================
rule summary_report:
    """Aggregate all scores into final leaderboard."""
    input:
        rules.run_evaluation.input
    output:
        csv=OUTPUT_DIR / "leaderboard.csv",
        html=OUTPUT_DIR / "leaderboard.html",
        json=OUTPUT_DIR / "summary.json"
    resources:
        cpus=1,
        mem_mb=4000,
        runtime=10
    shell:
        """
        python scripts/aggregate.py \
            --eval-dir {OUTPUT_DIR}/eval \
            --output-csv {output.csv} \
            --output-html {output.html} \
            --output-json {output.json} \
            --config config.yaml
        """

# =============================================================================
# Utility Rules
# =============================================================================
rule build_containers:
    """Build all Apptainer containers from Docker sources."""
    output:
        af3=CONTAINERS / "alphafold3.sif",
        protenix=CONTAINERS / "protenix.sif",
        colabfold=CONTAINERS / "colabfold.sif",
        mmseqs=CONTAINERS / "mmseqs2.sif",
        hhsuite=CONTAINERS / "hhsuite.sif"
    shell:
        """
        mkdir -p {CONTAINERS}
        
        # Build core inference containers (Docker build -> Apptainer)
        ./scripts/build_docker_containers.sh alphafold3
        ./scripts/build_docker_containers.sh protenix
        ./scripts/build_docker_containers.sh colabfold
        
        # Build tool containers (Direct pull from Docker registry)
        ./scripts/build_containers.sh mmseqs2
        ./scripts/build_containers.sh hhsuite
        """

rule clean_outputs:
    """Remove only generated predictions and reports."""
    shell:
        """
        rm -rf {OUTPUT_DIR}/pred {OUTPUT_DIR}/eval
        rm -f {OUTPUT_DIR}/summary_report.*
        """

rule clean_msa:
    """Remove MSA cache (expensive to regenerate)."""
    shell:
        """
        echo "Removing MSA cache in {CACHE_DIR}/sequences..."
        rm -rf {CACHE_DIR}/sequences
        """

rule clean_registry:
    """Remove registry database and generated target lists."""
    shell:
        """
        rm -rf {TARGET_LISTS}
        rm -f {REGISTRY_DB}
        rm -f {OUTPUT_DIR}/.sync_done
        """

rule clean_all:
    """[DANGEROUS] Remove all outputs and caches, but PRESERVE containers and manual test data."""
    shell:
        """
        rm -rf {OUTPUT_DIR} {CACHE_DIR} data/target_lists data/registry.sqlite
        """

# =============================================================================
# Default Target
# =============================================================================
rule all:
    input:
        rules.summary_report.output
    default_target: True
