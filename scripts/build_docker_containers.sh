#!/bin/bash
# =============================================================================
# Build Docker-based containers and convert to Apptainer
# =============================================================================
# Use this for containers that require Docker build (AlphaFold3, ColabFold, Protenix)
#
# Usage:
#   ./scripts/build_docker_containers.sh [alphafold3|colabfold|protenix|all]
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DOCKERFILE_DIR="${PROJECT_DIR}/containers/dockerfiles"
CONTAINER_DIR="${PROJECT_DIR}/containers"
BUILD_DIR="/tmp/psp-container-build"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'

log() { echo -e "${BLUE}[INFO]${NC} $1"; }
success() { echo -e "${GREEN}[OK]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; }

build_alphafold3() {
    log "Building AlphaFold3..."
    
    # Clone repo (AF3 Dockerfile needs to run from within repo)
    rm -rf "${BUILD_DIR}/alphafold3"
    git clone --depth 1 https://github.com/google-deepmind/alphafold3.git "${BUILD_DIR}/alphafold3"
    
    # Build Docker image using official Dockerfile
    cd "${BUILD_DIR}/alphafold3"
    docker build --network=host -t alphafold3-local -f docker/Dockerfile .
    
    # Convert to Apptainer
    log "Converting to Apptainer..."
    apptainer build "${CONTAINER_DIR}/alphafold3.sif" docker-daemon://alphafold3-local:latest
    
    success "Built: ${CONTAINER_DIR}/alphafold3.sif"
    ls -lh "${CONTAINER_DIR}/alphafold3.sif"
}

build_colabfold() {
    log "Building ColabFold..."
    
    local dockerfile="${DOCKERFILE_DIR}/Dockerfile.colabfold"
    
    docker build -t colabfold-local -f "$dockerfile" "$PROJECT_DIR"
    
    log "Converting to Apptainer..."
    apptainer build "${CONTAINER_DIR}/colabfold.sif" docker-daemon://colabfold-local:latest
    
    success "Built: ${CONTAINER_DIR}/colabfold.sif"
    ls -lh "${CONTAINER_DIR}/colabfold.sif"
}

build_protenix() {
    log "Building Protenix..."
    
    local dockerfile="${DOCKERFILE_DIR}/Dockerfile.protenix"
    
    docker build --network=host -t protenix-local -f "$dockerfile" "$PROJECT_DIR"
    
    log "Converting to Apptainer..."
    apptainer build "${CONTAINER_DIR}/protenix.sif" docker-daemon://protenix-local:latest
    
    success "Built: ${CONTAINER_DIR}/protenix.sif"
    ls -lh "${CONTAINER_DIR}/protenix.sif"
}

# Cleanup function
cleanup() {
    rm -rf "$BUILD_DIR"
}

# Main
mkdir -p "$CONTAINER_DIR" "$BUILD_DIR"

case "${1:-help}" in
    alphafold3)
        build_alphafold3
        ;;
    colabfold)
        build_colabfold
        ;;
    protenix)
        build_protenix
        ;;
    all)
        build_alphafold3 || error "AlphaFold3 failed"
        build_colabfold || error "ColabFold failed"
        build_protenix || error "Protenix failed"
        ;;
    *)
        echo "Usage: $0 [alphafold3|colabfold|protenix|all]"
        echo ""
        echo "Builds inference containers via Docker, then converts to Apptainer."
        echo "Requires: Docker with GPU support, Apptainer"
        exit 1
        ;;
esac

success "Done!"
