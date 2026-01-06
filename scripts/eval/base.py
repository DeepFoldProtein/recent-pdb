#!/usr/bin/env python3
"""
Evaluation Scorer Plugin Base
=============================
Abstract interface for structure comparison scorers.

To add a new scorer:
1. Create a new file in scripts/eval/ (e.g., my_scorer.py)
2. Inherit from ScorerPlugin
3. Implement the score() method
4. Add the scorer name to config.yaml enabled_scorers
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Optional
import json
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


@dataclass
class EvalResult:
    """
    Result from a single evaluation comparison.

    Attributes:
        scorer_name: Name of the scorer that produced this result
        target_id: PDB ID of the target
        model_name: Name of the prediction model
        score: Primary score value (higher is better)
        score_details: Additional per-residue or per-chain scores
        metadata: Any additional information
    """

    scorer_name: str
    target_id: str
    model_name: str
    score: float
    score_details: Optional[dict[str, Any]] = None
    metadata: Optional[dict[str, Any]] = None

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {k: v for k, v in asdict(self).items() if v is not None}

    def save(self, output_path: Path):
        """Save result to JSON file."""
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)
        logger.info(f"Saved result to {output_path}")


class ScorerPlugin(ABC):
    """
    Abstract base class for evaluation scorer plugins.

    Each scorer computes a specific metric comparing predicted
    structures to experimental ground truth.
    """

    # Override in subclass
    name: str = "base"
    description: str = "Base scorer plugin"

    # Score interpretation
    higher_is_better: bool = True
    score_range: tuple[float, float] = (0.0, 1.0)

    def __init__(self, config: dict):
        """
        Initialize scorer with configuration.

        Args:
            config: Full pipeline configuration dictionary
        """
        self.config = config

    @abstractmethod
    def score(
        self,
        prediction_path: Path,
        ground_truth_path: Path,
        target_id: str,
        model_name: str,
    ) -> EvalResult:
        """
        Compute evaluation metric between prediction and ground truth.

        Args:
            prediction_path: Path to predicted structure (CIF/PDB)
            ground_truth_path: Path to experimental structure (CIF/PDB)
            target_id: PDB ID of the target
            model_name: Name of the prediction model

        Returns:
            EvalResult with computed score and details
        """
        pass

    def validate_structures(
        self,
        prediction_path: Path,
        ground_truth_path: Path,
    ) -> bool:
        """Validate that both structure files exist and are readable."""
        if not prediction_path.exists():
            logger.error(f"Prediction file not found: {prediction_path}")
            return False
        if not ground_truth_path.exists():
            logger.error(f"Ground truth file not found: {ground_truth_path}")
            return False
        if prediction_path.stat().st_size == 0:
            logger.error(f"Prediction file is empty: {prediction_path}")
            return False
        if ground_truth_path.stat().st_size == 0:
            logger.error(f"Ground truth file is empty: {ground_truth_path}")
            return False
        return True


def load_structure_gemmi(path: Path):
    """
    Load structure using gemmi library.

    Args:
        path: Path to CIF or PDB file

    Returns:
        gemmi.Structure object
    """
    try:
        import gemmi
    except ImportError:
        raise ImportError("gemmi is required for structure loading")

    if path.suffix.lower() in (".cif", ".mmcif"):
        return gemmi.read_structure(str(path))
    else:
        return gemmi.read_pdb(str(path))


def get_ca_coordinates(structure) -> list[tuple[str, int, tuple[float, float, float]]]:
    """
    Extract CA atom coordinates from structure.

    Returns:
        List of (chain_id, residue_number, (x, y, z)) tuples
    """
    coords = []

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.name == "CA":
                        coords.append(
                            (
                                chain.name,
                                residue.seqid.num,
                                (atom.pos.x, atom.pos.y, atom.pos.z),
                            )
                        )
                        break  # Only one CA per residue
        break  # Only first model

    return coords
