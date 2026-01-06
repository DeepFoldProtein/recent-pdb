#!/bin/bash
# Setup resources for inference (one-time setup)

# Source environment variables if available
if [ -f "setup_env.sh" ]; then
    source setup_env.sh
else
    echo "Error: setup_env.sh not found. Please copy setup_env.sh.template to setup_env.sh and configure it."
    exit 1
fi

# 1. Protenix Resources (CCD, etc.)
echo "Setting up Protenix resources..."
mkdir -p "$PROTENIX_CACHE_DIR"

# Check if cache is already populated (simple check for ccd_cache)
if [ -d "$PROTENIX_CACHE_DIR/ccd_cache" ]; then
    echo "Protenix cache appears to exist in $PROTENIX_CACHE_DIR. Skipping download."
else
    # Run download command with PYTHONPATH setup
    # Note: runner is a top-level module in /opt/protenix, not sub-package of protenix
    # We manually construct the config since protenix.config doesn't expose get_config directly in this version
    # We explicitly import data_configs and inference_configs to populate required keys
    apptainer exec --nv \
      --bind "$PROTENIX_CACHE_DIR:/opt/protenix/release_data" \
      --env PYTHONPATH=/opt/protenix \
      "$PROTENIX_SIF" \
      python3 -c "
import sys; sys.path.append('/opt/protenix')
from ml_collections import ConfigDict
from configs.configs_base import configs
from configs.configs_data import data_configs
from configs.configs_inference import inference_configs
c = ConfigDict(configs)
c.data = ConfigDict(data_configs)
c.update(ConfigDict(inference_configs))
from runner.inference import download_infercence_cache
download_infercence_cache(c)
"
fi

# 2. ColabFold Weights
echo "Setting up ColabFold resources..."
mkdir -p "$COLABFOLD_CACHE_DIR"

# Check marker file (colabfold touches a file upon completion)
if [ -f "$COLABFOLD_CACHE_DIR/colabfold/params/download_complexes_multimer_v3_finished.txt" ]; then
     echo "ColabFold weights appear to exist. Skipping download."
else
    apptainer exec --nv \
      --bind "$COLABFOLD_CACHE_DIR:/cache" \
      "$COLABFOLD_SIF" \
      python3 -m colabfold.download
fi

echo "Resource setup complete."
