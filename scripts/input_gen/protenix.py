#!/usr/bin/env python3
"""
Protenix Input Generator
========================
Generates pickle (PKL) input format for Protenix inference.

Protenix uses a pickle-serialized dictionary with:
- sequences: Chain sequences
- features: Pre-computed MSA features
"""

import argparse
from pathlib import Path

import json

from .base import (
    ModelInputGenerator,
    TargetFeatures,
    load_target_features,
    logger,
)


class ProtenixInputGenerator(ModelInputGenerator):
    """Generate Protenix PKL input format."""

    name = "Protenix"
    input_format = "json"

    def generate(
        self,
        target_features: TargetFeatures,
        output_dir: Path,
    ) -> Path:
        """Generate Protenix input JSON."""
        output_dir.mkdir(parents=True, exist_ok=True)
        # Note: Protenix expects a list of jobs, even for one job
        output_path = output_dir / f"{target_features.target_id}.json"

        # Build sequence list
        sequences = []
        for chain in target_features.chains:
            # Basic protein chain definition
            protein_chain = {
                "proteinChain": {
                    "sequence": chain.sequence.replace("-", "")
                    .replace(".", "")
                    .upper(),
                    "count": 1,
                    # We leave 'msa' empty or unspecified to let Protenix handle it
                    # or if we want to provide precomputed, we need to restructure them.
                    # For this test, we rely on Protenix's internal search or if it supports
                    # external MSAs differently (e.g. via CLI flags).
                    # The user's goal implies we might want to use our MSAs, but Protenix's
                    # JSON format requires specific directory structure (pairing.a3m).
                    # For now, we generate the valid structure input.
                }
            }
            sequences.append(protein_chain)

        job_entry = {
            "name": target_features.target_id,
            "sequences": sequences,
        }

        # Protenix input is a list of job dictionaries
        protenix_input = [job_entry]

        # Write JSON
        with open(output_path, "w") as f:
            json.dump(protenix_input, f, indent=2)

        logger.info(f"Generated Protenix input: {output_path}")
        return output_path

    def _encode_msa(self, *args, **kwargs):
        # Deprecated for JSON format
        pass


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
    parser = argparse.ArgumentParser(description="Generate Protenix input pickle")
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

    generator = ProtenixInputGenerator(config)
    output_path = generator.generate(target_features, args.output_dir / args.target_id)

    if generator.validate_output(output_path):
        logger.info("Input generation successful")
        return 0
    else:
        logger.error("Input generation failed")
        return 1


if __name__ == "__main__":
    exit(main())
