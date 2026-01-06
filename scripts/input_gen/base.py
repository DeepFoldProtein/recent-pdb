#!/usr/bin/env python3
"""
Model Input Generator Base Class
================================
Abstract interface for generating model-specific input features from shared MSA cache.

To add a new model:
1. Create a new file in scripts/input_gen/ (e.g., mymodel.py)
2. Inherit from ModelInputGenerator
3. Implement the generate() method
4. Add the model name to config.yaml enabled_models
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import json
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


@dataclass
class ChainFeatures:
    """Features for a single chain."""

    chain_id: str
    sequence: str
    seq_hash: str
    msa_paths: list[Path] = None
    template_paths: Optional[list[Path]] = None


@dataclass
class TargetFeatures:
    """Aggregated features for a prediction target."""

    target_id: str
    chains: list[ChainFeatures]

    @property
    def num_chains(self) -> int:
        return len(self.chains)

    @property
    def total_length(self) -> int:
        return sum(len(c.sequence) for c in self.chains)


class ModelInputGenerator(ABC):
    """
    Abstract base class for model-specific input generation.

    Each model (AlphaFold3, Protenix, RoseTTAFold2NA, etc.) has different
    input format requirements. This class provides a common interface.
    """

    # Override in subclass
    name: str = "base"
    input_format: str = "unknown"  # json, pkl, npz, etc.

    def __init__(self, config: dict):
        """
        Initialize generator with configuration.

        Args:
            config: Full pipeline configuration dictionary
        """
        self.config = config
        self.model_config = config.get("models", {}).get(self.name, {})

    @abstractmethod
    def generate(
        self,
        target_features: TargetFeatures,
        output_dir: Path,
    ) -> Path:
        """
        Generate model-specific input from target features.

        Args:
            target_features: Aggregated features including MSA paths
            output_dir: Directory to write output files

        Returns:
            Path to the main generated input file
        """
        pass

    def load_msa(self, msa_path: Path) -> list[tuple[str, str]]:
        """
        Load MSA from A3M file.

        Args:
            msa_path: Path to A3M file

        Returns:
            List of (name, sequence) tuples
        """
        sequences = []

        if not msa_path.exists():
            logger.warning(f"MSA file not found: {msa_path}")
            return sequences

        with open(msa_path) as f:
            current_name = None
            current_seq = []

            for line in f:
                if line.startswith(">"):
                    if current_name:
                        sequences.append((current_name, "".join(current_seq)))
                    current_name = line.strip()[1:]
                    current_seq = []
                else:
                    current_seq.append(line.strip())

            if current_name:
                sequences.append((current_name, "".join(current_seq)))

        return sequences

    def validate_output(self, output_path: Path) -> bool:
        """Validate that output file was generated correctly."""
        if not output_path.exists():
            return False
        if output_path.stat().st_size == 0:
            return False
        return True


def load_target_features(
    target_id: str,
    cache_dir: Path,
    target_lists_dir: Path,
    config: dict = None,
    featurizer: str = None,
) -> Optional[TargetFeatures]:
    """
    Load target features from registry and cache.

    Args:
        target_id: PDB ID of target
        cache_dir: Shared MSA cache directory
        target_lists_dir: Directory containing target details
        config: Full pipeline config
        featurizer: Featurizer name to filter MSA sources

    Returns:
        TargetFeatures object or None if not found
    """
    details_file = target_lists_dir / "target_details.json"

    if not details_file.exists():
        logger.error(f"Target details not found: {details_file}")
        return None

    with open(details_file) as f:
        all_targets = json.load(f)

    # Find target
    target_info = None
    for t in all_targets:
        if t["pdb_id"].upper() == target_id.upper():
            target_info = t
            break

    if target_info is None:
        logger.error(f"Target not found in registry: {target_id}")
        return None

    # Get featurizer config if specified
    featurizer_sources = []
    if featurizer and config:
        featurizer_cfg = config.get("featurizers", {}).get(featurizer, {})
        featurizer_sources = featurizer_cfg.get("sources", [])

    # Build chain features
    chains = []
    for chain_info in target_info["chains"]:
        seq_hash = chain_info["seq_hash"]

        # Collect MSA paths based on featurizer sources
        msa_paths = []

        if featurizer_sources:
            # Use only sources specified in featurizer config
            for src in featurizer_sources:
                aligner = src["aligner"]
                db_id = src["database"]

                # Get tool version from config
                tool_cfg = config.get("msa", {}).get("tools", {}).get(aligner, {})
                version = tool_cfg.get("version", "")
                tool_dir_name = f"{aligner}_{version}" if version else aligner

                # Get params hash
                import hashlib

                params = tool_cfg.get("params", {})
                if params is None:
                    params = {}
                params_str = json.dumps(params, sort_keys=True)
                phash = hashlib.md5(params_str.encode()).hexdigest()[:8]

                # Build path
                msa_file = cache_dir / tool_dir_name / db_id / phash / f"{seq_hash}.a3m"
                if msa_file.exists():
                    msa_paths.append(msa_file)
                else:
                    logger.warning(f"MSA not found: {msa_file}")
        else:
            # Fallback: collect all available MSAs
            for tool_prefix in ["mmseqs2", "hhblits"]:
                for tool_dir in cache_dir.glob(f"{tool_prefix}*"):
                    if not tool_dir.is_dir():
                        continue
                    for db_dir in tool_dir.iterdir():
                        if not db_dir.is_dir():
                            continue
                        for param_dir in db_dir.iterdir():
                            if not param_dir.is_dir():
                                continue
                            msa_file = param_dir / f"{seq_hash}.a3m"
                            if msa_file.exists():
                                msa_paths.append(msa_file)

        if not msa_paths:
            # Fallback/Legacy: check merged
            merged = cache_dir / seq_hash / "msa" / "merged.a3m"
            if merged.exists():
                msa_paths.append(merged)

        chains.append(
            ChainFeatures(
                chain_id=chain_info["chain_id"],
                sequence=chain_info["sequence"],
                seq_hash=seq_hash,
                msa_paths=msa_paths,
            )
        )

    return TargetFeatures(
        target_id=target_id,
        chains=chains,
    )
