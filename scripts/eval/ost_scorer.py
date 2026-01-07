#!/usr/bin/env python3
"""
OST (OpenStructure) Scorer Plugin
=================================
Computes lDDT (local Distance Difference Test) using OpenStructure.

lDDT is a per-residue score measuring local structural similarity,
commonly used in CASP and protein structure prediction benchmarks.
"""

import argparse
import json
import subprocess
import tempfile
from pathlib import Path


from .base import ScorerPlugin, EvalResult, logger

# gemmi is used for SEQRES parsing in the host environment
try:
    import gemmi

    HAS_GEMMI = True
except ImportError:
    HAS_GEMMI = False


class OSTScorer(ScorerPlugin):
    """OpenStructure-based lDDT scorer using Apptainer."""

    name = "ost"
    description = "lDDT score using OpenStructure (via Apptainer)"
    higher_is_better = True
    score_range = (0.0, 1.0)

    def score(
        self,
        prediction_path: Path,
        ground_truth_path: Path,
        target_id: str,
        model_name: str,
    ) -> EvalResult:
        """Compute lDDT between prediction and ground truth."""
        if not self.validate_structures(prediction_path, ground_truth_path):
            return EvalResult(
                scorer_name=self.name,
                target_id=target_id,
                model_name=model_name,
                score=0.0,
                metadata={"error": "Structure validation failed"},
            )

        container_path = self.config.get("containers", {}).get("openstructure")
        if not container_path:
            # Try default path
            container_path = "containers/openstructure.sif"

        if Path(container_path).exists():
            return self._score_with_container(
                prediction_path,
                ground_truth_path,
                target_id,
                model_name,
                container_path,
            )
        else:
            logger.warning(
                f"OST container not found at {container_path}. Using fallback."
            )
            return self._score_fallback(
                prediction_path, ground_truth_path, target_id, model_name
            )

    def _parse_seqres_from_cif(
        self, cif_path: Path
    ) -> tuple[dict[str, str], dict[str, str]]:
        """
        Parse SEQRES from mmCIF file using gemmi.
        Uses _entity_poly_seq for sequence and _struct_asym for chain mapping.

        Returns:
            - mapping: dict mapping entity_id to one_letter_code
            - chain_mapping: dict mapping chain_id (asym_id) to entity_id
        """
        if not HAS_GEMMI:
            raise ImportError("gemmi is required for SEQRES parsing")

        doc = gemmi.cif.read_file(str(cif_path))
        block = doc.sole_block()

        # 1. Map entity_id to sequence using _entity_poly_seq
        # This is more robust than _entity_poly.pdbx_seq_one_letter_code
        entity_seqs = {}
        poly_seq = block.get_mmcif_category("_entity_poly_seq.")
        if "entity_id" in poly_seq and "mon_id" in poly_seq:
            # Group by entity_id
            from collections import defaultdict

            eid_to_res = defaultdict(list)
            for eid, mon in zip(poly_seq["entity_id"], poly_seq["mon_id"]):
                # Convert 3-letter to 1-letter
                one_letter = gemmi.one_letter_code([mon])
                # If unknown, use 'X'
                if not one_letter or one_letter == " ":
                    one_letter = "X"
                eid_to_res[eid].append(one_letter)

            for eid, res_list in eid_to_res.items():
                entity_seqs[eid] = "".join(res_list)

        # Fallback to _entity_poly if _entity_poly_seq is missing
        if not entity_seqs:
            poly = block.get_mmcif_category("_entity_poly.")
            if "entity_id" in poly and "pdbx_seq_one_letter_code" in poly:
                for eid, seq in zip(
                    poly["entity_id"], poly["pdbx_seq_one_letter_code"]
                ):
                    clean_seq = "".join(seq.split())
                    entity_seqs[eid] = clean_seq

        # 2. Map chain_id (asym_id) to entity_id
        struct_asym = block.get_mmcif_category("_struct_asym.")
        chain_to_entity = {}
        if "id" in struct_asym and "entity_id" in struct_asym:
            for asym_id, eid in zip(struct_asym["id"], struct_asym["entity_id"]):
                if eid in entity_seqs:
                    chain_to_entity[asym_id] = eid

        return entity_seqs, chain_to_entity

    def _score_with_container(
        self,
        prediction_path: Path,
        ground_truth_path: Path,
        target_id: str,
        model_name: str,
        container_path: str,
    ) -> EvalResult:
        """Score using OpenStructure 'compare-structures' in Apptainer."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            scores_json = tmp_path / "scores.json"
            seqres_fasta = tmp_path / "seqres.fasta"

            # Parse SEQRES if reference is CIF
            mapping_args = []
            if ground_truth_path.suffix.lower() in (".cif", ".mmcif"):
                try:
                    entity_seqs, chain_to_entity = self._parse_seqres_from_cif(
                        ground_truth_path
                    )
                    if entity_seqs:
                        with open(seqres_fasta, "w") as f:
                            for eid, seq in entity_seqs.items():
                                f.write(f">{eid}\n{seq}\n")

                        mapping_args.append("--seqres")
                        mapping_args.append(str(seqres_fasta))
                        mapping_args.append("--trg-seqres-mapping")
                        for asym_id, eid in chain_to_entity.items():
                            mapping_args.append(f"{asym_id}:{eid}")
                        mapping_args.append("--residue-number-alignment")
                except Exception as e:
                    logger.warning(
                        f"Failed to parse SEQRES from {ground_truth_path}: {e}"
                    )

            # Prepare command
            # Bind necessary directories: current, reference, and model
            binds = [f"{Path().absolute()}:/work"]

            # Add bind for reference if it's not under current dir
            ref_abs = ground_truth_path.absolute()
            pred_abs = prediction_path.absolute()

            # Helper to find common mount points (simplistic for /store)
            if str(ref_abs).startswith("/store"):
                binds.append("/store:/store")
            elif str(ref_abs).startswith("/home"):
                binds.append("/home:/home")

            # Ensure model path is also bound if different
            if str(pred_abs).startswith("/store") and "/store:/store" not in binds:
                binds.append("/store:/store")

            cmd = [
                "apptainer",
                "exec",
            ]
            for b in set(binds):
                cmd.extend(["--bind", b])

            cmd.extend(
                [
                    container_path,
                    "ost",
                    "compare-structures",
                    "--model",
                    str(pred_abs),
                    "--reference",
                    str(ref_abs),
                    "--output",
                    str(scores_json),
                    "--lddt",
                    "--local-lddt",
                    "--qs-score",
                    "--dockq",
                    "--ics",
                    "--ips",
                    "--tm-score",
                ]
            )
            cmd += mapping_args

            try:
                subprocess.run(
                    cmd,
                    check=True,
                    capture_output=True,
                    text=True,
                )

                if not scores_json.exists():
                    raise FileNotFoundError(f"OST output {scores_json} not found")

                with open(scores_json) as f:
                    data = json.load(f)

                # Extract global lDDT
                # OST JSON format: {"lddt": 0.85, "local_lddt": {"A": [0.8, ...]}, ...}
                global_lddt = data.get("lddt", 0.0)
                qs_global = data.get("qs_global", 0.0)
                dockq = data.get("dockq", 0.0)
                ics = data.get("ics", 0.0)
                ips = data.get("ips", 0.0)
                tm_score = data.get("tm_score", 0.0)

                return EvalResult(
                    scorer_name=self.name,
                    target_id=target_id,
                    model_name=model_name,
                    score=global_lddt,
                    score_details=data,
                    metadata={
                        "container": container_path,
                        "used_seqres": len(mapping_args) > 0,
                        "qs_score": qs_global,
                        "dockq": dockq,
                        "ics": ics,
                        "ips": ips,
                        "tm_score": tm_score,
                    },
                )

            except subprocess.CalledProcessError as e:
                logger.error(f"OST container failed: {e.stderr}")
                return EvalResult(
                    scorer_name=self.name,
                    target_id=target_id,
                    model_name=model_name,
                    score=0.0,
                    metadata={"error": e.stderr},
                )
            except Exception as e:
                logger.error(f"OST scoring failed: {e}")
                return EvalResult(
                    scorer_name=self.name,
                    target_id=target_id,
                    model_name=model_name,
                    score=0.0,
                    metadata={"error": str(e)},
                )

    def _score_fallback(
        self,
        prediction_path: Path,
        ground_truth_path: Path,
        target_id: str,
        model_name: str,
    ) -> EvalResult:
        """Simplified lDDT calculation without OST."""
        try:
            from .base import load_structure_gemmi, get_ca_coordinates
            import numpy as np

            # Load structures
            pred_struct = load_structure_gemmi(prediction_path)
            ref_struct = load_structure_gemmi(ground_truth_path)

            pred_coords = get_ca_coordinates(pred_struct)
            ref_coords = get_ca_coordinates(ref_struct)

            # Match residues by chain and number
            pred_dict = {(c, n): xyz for c, n, xyz in pred_coords}
            ref_dict = {(c, n): xyz for c, n, xyz in ref_coords}

            common_keys = set(pred_dict.keys()) & set(ref_dict.keys())

            if len(common_keys) < 10:
                return EvalResult(
                    scorer_name=self.name,
                    target_id=target_id,
                    model_name=model_name,
                    score=0.0,
                    metadata={"error": "Too few common residues"},
                )

            # Compute pairwise distances
            keys = sorted(common_keys)
            pred_xyz = np.array([pred_dict[k] for k in keys])
            ref_xyz = np.array([ref_dict[k] for k in keys])

            n = len(keys)
            cutoffs = [0.5, 1.0, 2.0, 4.0]

            # Compute lDDT
            preserved = 0
            total = 0

            for i in range(n):
                for j in range(i + 1, n):
                    ref_dist = np.linalg.norm(ref_xyz[i] - ref_xyz[j])
                    if ref_dist > 15.0:  # Only consider local distances
                        continue

                    pred_dist = np.linalg.norm(pred_xyz[i] - pred_xyz[j])
                    diff = abs(ref_dist - pred_dist)

                    for cutoff in cutoffs:
                        total += 1
                        if diff < cutoff:
                            preserved += 1

            lddt = preserved / max(total, 1)

            return EvalResult(
                scorer_name=self.name,
                target_id=target_id,
                model_name=model_name,
                score=lddt,
                score_details={
                    "num_residues": n,
                    "num_pairs": total // len(cutoffs),
                },
                metadata={
                    "cutoffs": cutoffs,
                    "library": "fallback",
                },
            )
        except Exception as e:
            logger.error(f"Fallback lDDT failed: {e}")
            return EvalResult(
                scorer_name=self.name,
                target_id=target_id,
                model_name=model_name,
                score=0.0,
                metadata={"error": str(e)},
            )


def _load_config(config_path):
    """Load config with local overrides."""
    import yaml
    from pathlib import Path
    config_path = Path(config_path)
    with open(config_path) as f:
        config = yaml.safe_load(f)
    local_path = config_path.parent / "config.local.yaml"
    if local_path.exists():
        with open(local_path) as f:
            local_config = yaml.safe_load(f) or {}
        for key, value in local_config.items():
            if key in config and isinstance(config[key], dict) and isinstance(value, dict):
                config[key].update(value)
            else:
                config[key] = value
    return config


def main():
    parser = argparse.ArgumentParser(description="Compute lDDT score")
    parser.add_argument("--prediction", type=Path, required=True)
    parser.add_argument("--ground-truth", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=Path("config.yaml"))

    parser.add_argument("--target-id", type=str)
    parser.add_argument("--model-name", type=str)

    args = parser.parse_args()

    # Extract target and model from paths or arguments
    target_id = args.target_id or args.prediction.stem
    model_name = args.model_name or (
        args.prediction.parent.parent.name
        if len(args.prediction.parents) > 1
        else "unknown"
    )

    config = _load_config(args.config)

    scorer = OSTScorer(config)
    result = scorer.score(
        args.prediction,
        args.ground_truth,
        target_id,
        model_name,
    )

    result.save(args.output)

    logger.info(f"lDDT = {result.score:.4f}")
    return 0 if result.score > 0 else 1


if __name__ == "__main__":
    exit(main())
