#!/usr/bin/env python3
"""
AlphaFold3 Input Generator
==========================
Generates JSON input format for AlphaFold3 inference.

AlphaFold3 expects a specific JSON schema with:
- sequences: List of chain sequences with entity information
- msas: MSA data per entity
- templates: Template structures (optional)
"""

import argparse

import json
from pathlib import Path


from .base import (
    ModelInputGenerator,
    TargetFeatures,
    load_target_features,
    logger,
)


class AlphaFold3InputGenerator(ModelInputGenerator):
    """Generate AlphaFold3 JSON input format."""

    name = "AlphaFold3"
    input_format = "json"

    def generate(
        self,
        target_features: TargetFeatures,
        output_dir: Path,
    ) -> Path:
        """Generate AlphaFold3 input JSON."""
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{target_features.target_id}.json"

        # Build input structure
        af3_input = {
            "name": target_features.target_id,
            "dialect": "alphafold3",
            "version": 1,
            "modelSeeds": [42],  # Default seed
            "sequences": [],
        }

        # Group chains by sequence hash (for homo-oligomers)
        hash_to_chains: dict[str, list] = {}
        for chain in target_features.chains:
            if chain.seq_hash not in hash_to_chains:
                hash_to_chains[chain.seq_hash] = []
            hash_to_chains[chain.seq_hash].append(chain)

        # Build sequence entries
        for seq_hash, chains in hash_to_chains.items():
            first_chain = chains[0]

            # Load and combine MSAs
            msa_data = []
            seen_seqs = set()
            templates = []

            # Add query sequence first
            if first_chain.sequence:
                seen_seqs.add(first_chain.sequence)
                msa_data.append({"sequence": first_chain.sequence})

            if first_chain.msa_paths:
                for path in first_chain.msa_paths:
                    # Check if this is a template database (pdb100)
                    if "pdb100" in str(path) and not templates:
                        templates = self._process_templates(path, first_chain.sequence)
                        continue

                    # Otherwise treat as MSA
                    if "pdb100" not in str(path):
                        msa_seqs = self.load_msa(path)
                        for _, seq in msa_seqs:
                            # Standardize A3M -> Aligned FASTA (remove insertions)
                            # Keep uppercase (matches) and dashes (deletions)
                            seq_clean = "".join(
                                c for c in seq if c.isupper() or c == "-"
                            )
                            if seq_clean not in seen_seqs:
                                seen_seqs.add(seq_clean)
                                msa_data.append({"sequence": seq_clean})

            # Limit MSA depth
            msa_data = msa_data[:2048]

            # Add MSA if available
            if msa_data:
                unpaired_msa = "\n".join(
                    f">{i}\n{m['sequence']}" for i, m in enumerate(msa_data)
                )
            else:
                unpaired_msa = f">{first_chain.sequence}\n{first_chain.sequence}"

            sequence_entry = {
                "protein": {
                    "id": [c.chain_id for c in chains],
                    "sequence": first_chain.sequence,
                    "unpairedMsa": unpaired_msa,
                    "pairedMsa": "",
                    "templates": templates,
                },
            }

            af3_input["sequences"].append(sequence_entry)

        # Write output
        with open(output_path, "w") as f:
            json.dump(af3_input, f, indent=2)

        logger.info(f"Generated AlphaFold3 input: {output_path}")
        return output_path

    def _process_templates(self, a3m_path: Path, query_sequence: str) -> list[dict]:
        """
        Process A3M file to find templates and generate mmCIFs.
        :param a3m_path: Path to A3M file
        :param query_sequence: Query sequence (used for validation/logging)
        :return: List of template objects
        """
        if not self.config.get("input_generation", {}).get("use_templates", True):
            logger.info("Template generation disabled in config")
            return []

        templates = []
        return templates


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
    parser = argparse.ArgumentParser(description="Generate AlphaFold3 input JSON")
    parser.add_argument("--target-id", required=True, help="Target PDB ID")
    parser.add_argument("--featurizer", required=True, help="Featurizer name")
    parser.add_argument(
        "--cache-dir", type=Path, required=True, help="MSA cache directory"
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True, help="Output directory"
    )
    parser.add_argument("--config", type=Path, default=Path("config.yaml"))

    args = parser.parse_args()

    config = _load_config(args.config)

    target_lists_dir = Path(config["paths"]["target_lists"])
    target_features = load_target_features(
        args.target_id,
        args.cache_dir,
        target_lists_dir,
        config,
        featurizer=args.featurizer,
    )

    if target_features is None:
        return 1

    generator = AlphaFold3InputGenerator(config)
    output_path = generator.generate(target_features, args.output_dir / args.target_id)

    if generator.validate_output(output_path):
        logger.info("Input generation successful")
        return 0
    else:
        logger.error("Input generation failed")
        return 1


if __name__ == "__main__":
    exit(main())
