#!/bin/bash
# =============================================================================
# Container Build Script for PSP Benchmark Pipeline
# =============================================================================
# Builds Apptainer containers from official Docker images.
#
# Usage:
#   ./scripts/build_containers.sh              # Build all containers
#   ./scripts/build_containers.sh alphafold3   # Build specific container
#   ./scripts/build_containers.sh --list       # List available containers
#
# References:
#   - AlphaFold3: https://github.com/google-deepmind/alphafold3
#   - Protenix:   https://github.com/bytedance/Protenix
#   - ColabFold:  https://github.com/sokrypton/ColabFold
#   - HH-suite:   https://github.com/soedinglab/hh-suite
#   - MMseqs2:    https://github.com/soedinglab/MMseqs2
#   - OpenStructure: https://openstructure.org/install
# =============================================================================

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
CONTAINER_DIR="${PROJECT_DIR}/containers"
LOG_DIR="${PROJECT_DIR}/logs"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Container definitions (directly pullable from registries)
# NOTE: For alphafold3, protenix, colabfold - use build_docker_containers.sh
declare -A CONTAINERS=(
    # MSA search containers
    ["mmseqs2"]="docker://ghcr.io/soedinglab/mmseqs2:latest"
    ["hhsuite"]="docker://soedinglab/hh-suite:latest"
    ["hmmer"]="docker://biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1"
    
    # Evaluation containers
    ["openstructure"]="docker://registry.scicore.unibas.ch/schwede/openstructure:latest"
)

# Container descriptions
declare -A DESCRIPTIONS=(
    ["mmseqs2"]="MMseqs2 sequence search (CPU/GPU)"
    ["hhsuite"]="HH-suite homology search (CPU)"
    ["hmmer"]="HMMER jackhmmer search (CPU)"
    ["openstructure"]="OpenStructure for lDDT scoring (CPU)"
)

# Helper functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

show_usage() {
    echo "Usage: $0 [OPTIONS] [CONTAINER...]"
    echo ""
    echo "Options:"
    echo "  --list, -l     List available containers"
    echo "  --force, -f    Force rebuild even if container exists"
    echo "  --fakeroot     Use fakeroot for building (no sudo)"
    echo "  --help, -h     Show this help message"
    echo ""
    echo "Containers:"
    for name in "${!CONTAINERS[@]}"; do
        printf "  %-15s %s\n" "$name" "${DESCRIPTIONS[$name]}"
    done
}

list_containers() {
    echo "Available containers:"
    echo "====================="
    for name in "${!CONTAINERS[@]}"; do
        local sif="${CONTAINER_DIR}/${name}.sif"
        if [[ -f "$sif" ]]; then
            local size=$(du -h "$sif" | cut -f1)
            echo -e "  ${GREEN}✓${NC} ${name} (${size})"
        else
            echo -e "  ${YELLOW}○${NC} ${name} (not built)"
        fi
    done
}

build_container() {
    local name="$1"
    local force="${2:-false}"
    local fakeroot="${3:-false}"
    
    if [[ ! -v "CONTAINERS[$name]" ]]; then
        log_error "Unknown container: $name"
        return 1
    fi
    
    local source="${CONTAINERS[$name]}"
    local output="${CONTAINER_DIR}/${name}.sif"
    local log_file="${LOG_DIR}/build_${name}.log"
    local def_file="${PROJECT_DIR}/containers/dockerfiles/${name}.def"
    
    # Check if already exists
    if [[ -f "$output" && "$force" != "true" ]]; then
        log_warn "$name.sif already exists. Use --force to rebuild."
        return 0
    fi
    
    # Prefer definition file if it exists
    if [[ -f "$def_file" ]]; then
        log_info "Building ${name} from definition file: ${def_file}"
        source="$def_file"
    else
        log_info "Building ${name} from ${source}..."
    fi
    
    # Build command
    local build_cmd="apptainer build --ignore-fakeroot-command"
    if [[ "$fakeroot" == "true" ]]; then
        build_cmd+=" --fakeroot"
    fi
    
    # Execute build
    mkdir -p "$CONTAINER_DIR" "$LOG_DIR"
    
    if $build_cmd "$output" "$source" 2>&1 | tee "$log_file"; then
        log_success "Built ${name}.sif successfully"
        ls -lh "$output"
    else
        log_error "Failed to build ${name}. Check ${log_file} for details."
        return 1
    fi
}

check_existing_containers() {
    log_info "Checking for existing containers in common locations..."
    
    # Known container locations
    local known_paths=(
        "/store/database/mmseqs/mmseqs.sif"
        "/shared/containers"
        "/opt/containers"
    )
    
    for path in "${known_paths[@]}"; do
        if [[ -f "$path" ]]; then
            log_info "Found: $path"
        elif [[ -d "$path" ]]; then
            log_info "Directory: $path"
            ls -la "$path"/*.sif 2>/dev/null || true
        fi
    done
}

# Link existing container instead of building
link_existing() {
    local name="$1"
    local existing_path="$2"
    
    if [[ ! -f "$existing_path" ]]; then
        log_error "Source file not found: $existing_path"
        return 1
    fi
    
    local target="${CONTAINER_DIR}/${name}.sif"
    mkdir -p "$CONTAINER_DIR"
    
    if [[ -f "$target" ]]; then
        log_warn "$target already exists"
        return 0
    fi
    
    ln -s "$existing_path" "$target"
    log_success "Linked $name.sif -> $existing_path"
}

# =============================================================================
# Main
# =============================================================================
main() {
    local force=false
    local fakeroot=false
    local containers_to_build=()
    
    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --list|-l)
                list_containers
                exit 0
                ;;
            --force|-f)
                force=true
                shift
                ;;
            --fakeroot)
                fakeroot=true
                shift
                ;;
            --help|-h)
                show_usage
                exit 0
                ;;
            --check)
                check_existing_containers
                exit 0
                ;;
            --link)
                # Usage: --link mmseqs2 /path/to/existing.sif
                shift
                link_existing "$1" "$2"
                exit 0
                ;;
            *)
                containers_to_build+=("$1")
                shift
                ;;
        esac
    done
    
    # Default to all containers if none specified
    if [[ ${#containers_to_build[@]} -eq 0 ]]; then
        containers_to_build=("${!CONTAINERS[@]}")
    fi
    
    # Create directories
    mkdir -p "$CONTAINER_DIR" "$LOG_DIR"
    
    # Build each container
    local failed=0
    for name in "${containers_to_build[@]}"; do
        if ! build_container "$name" "$force" "$fakeroot"; then
            ((failed++)) || true
        fi
    done
    
    echo ""
    log_info "Build summary:"
    list_containers
    
    if [[ $failed -gt 0 ]]; then
        log_error "$failed container(s) failed to build"
        exit 1
    fi
}

# Run if executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
