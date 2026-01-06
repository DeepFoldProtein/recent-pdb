#!/bin/bash
# =============================================================================
# Resource Management Script for PSP Benchmark
# =============================================================================
# Manages persistent resources (model weights, CCD files) in resources/
#
# Usage:
#   ./scripts/manage_resources.sh --status    # Check resource status
#   ./scripts/manage_resources.sh --clean     # [DANGEROUS] Remove all persistent resources
# =============================================================================

set -euo pipefail

# Source environment
if [ -f "setup_env.sh" ]; then
    source setup_env.sh
else
    echo "Error: setup_env.sh not found."
    exit 1
fi

# Configuration
RESOURCES_DIR="./resources"

show_status() {
    echo "Resource Status:"
    echo "================"
    
    if [ -d "$PROTENIX_CACHE_DIR" ]; then
        size=$(du -sh "$PROTENIX_CACHE_DIR" | cut -f1)
        echo -e "Protenix Cache: [EXIST] ($size) -> $PROTENIX_CACHE_DIR"
    else
        echo -e "Protenix Cache: [MISSING] -> $PROTENIX_CACHE_DIR"
    fi

    if [ -d "$COLABFOLD_CACHE_DIR" ]; then
        size=$(du -sh "$COLABFOLD_CACHE_DIR" | cut -f1)
        echo -e "ColabFold Cache: [EXIST] ($size) -> $COLABFOLD_CACHE_DIR"
    else
        echo -e "ColabFold Cache: [MISSING] -> $COLABFOLD_CACHE_DIR"
    fi
}

clean_resources() {
    echo "WARNING: This will delete all persistent model weights and CCD files in $RESOURCES_DIR."
    read -p "Are you sure you want to continue? (y/N) " confirm
    if [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]]; then
        echo "Cleaning resources..."
        rm -rf "$RESOURCES_DIR"/*
        echo "Resources cleaned."
    else
        echo "Cleanup cancelled."
    fi
}

case "${1:-status}" in
    --status)
        show_status
        ;;
    --clean)
        clean_resources
        ;;
    *)
        echo "Usage: $0 [--status|--clean]"
        exit 1
        ;;
esac
