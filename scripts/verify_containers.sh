#!/bin/bash
# =============================================================================
# Container Verification Script
# =============================================================================
# Checks if the built Apptainer containers can execute basic commands.
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
CONTAINER_DIR="${PROJECT_DIR}/containers"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

log_success() { echo -e "${GREEN}[OK]${NC} $1"; }
log_fail() { echo -e "${RED}[FAIL]${NC} $1"; }

verify_container() {
    local name="$1"
    local cmd="$2"
    local expected="$3"
    local sif="${CONTAINER_DIR}/${name}.sif"

    echo -n "Verifying ${name}..."
    
    if [[ ! -f "$sif" ]]; then
        echo -e " ${RED}Not found${NC}"
        return
    fi

    if apptainer exec "$sif" $cmd 2>&1 | grep -i -q "$expected"; then
        echo -e " ${GREEN}OK${NC}"
    else
        echo -e " ${RED}Failed${NC}"
        echo "Command: apptainer exec $sif $cmd"
        apptainer exec "$sif" $cmd 2>&1 | head -n 5
    fi
}

echo "Checking containers in ${CONTAINER_DIR}..."
echo "----------------------------------------"

# MSA Tools
verify_container "mmseqs2" "mmseqs -h" "MMseqs2"
verify_container "hhsuite" "hhblits -h" "HHblits"
verify_container "hmmer" "jackhmmer -h" "jackhmmer"

# Evaluation Tools
verify_container "openstructure" "ost -h" "usage"

# Inference Tools (if built)
verify_container "alphafold3" "python /app/alphafold/run_alphafold.py --help" "AlphaFold 3 structure prediction script"
verify_container "colabfold" "colabfold_batch --help" "usage: colabfold_batch"
verify_container "protenix" "protenix --help" "usage: protenix"

echo "----------------------------------------"
