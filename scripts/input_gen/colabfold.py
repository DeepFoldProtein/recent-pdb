#!/usr/bin/env python3
"""
ColabFold Input Generator
=========================
Generates FASTA or A3M input format for ColabFold inference.

ColabFold accepts:
- FASTA files for automatic MSA generation
- A3M files for pre-computed MSAs
- Can output AF3-compatible JSON with --af3-json flag
"""

import argparse
from pathlib import Path


from .base import (
    ModelInputGenerator,
    TargetFeatures,
    load_target_features,
    logger,
)


class ColabFoldInputGenerator(ModelInputGenerator):
    """Generate ColabFold FASTA/A3M input format."""

    name = "ColabFold"
    input_format = "fasta"

    def generate(
        self,
        target_features: TargetFeatures,
        output_dir: Path,
    ) -> Path:
        """Generate ColabFold input FASTA with MSA if available."""
        output_dir.mkdir(parents=True, exist_ok=True)

        # Check if we have pre-computed MSAs
        has_msa = any(c.msa_paths for c in target_features.chains)

        if has_msa:
            # Generate A3M input with pre-computed MSAs
            return self._generate_a3m(target_features, output_dir)
        else:
            # Generate FASTA for ColabFold to compute MSAs
            return self._generate_fasta(target_features, output_dir)

    def _generate_fasta(
        self,
        target_features: TargetFeatures,
        output_dir: Path,
    ) -> Path:
        """Generate FASTA file for ColabFold."""
        output_path = output_dir / f"{target_features.target_id}.fasta"

        with open(output_path, "w") as f:
            if target_features.num_chains == 1:
                # Single chain - simple format
                chain = target_features.chains[0]
                f.write(f">{target_features.target_id}\n{chain.sequence}\n")
            else:
                # Multi-chain complex - use ColabFold complex format
                # Format: >name|chain1:chain2
                chain_names = [f"chain{c.chain_id}" for c in target_features.chains]
                sequences = [c.sequence for c in target_features.chains]

                header = f">{target_features.target_id}|{':'.join(chain_names)}"
                seq_line = ":".join(sequences)

                f.write(f"{header}\n{seq_line}\n")

        logger.info(f"Generated ColabFold FASTA: {output_path}")
        return output_path

    def _generate_a3m(
        self,
        target_features: TargetFeatures,
        output_dir: Path,
    ) -> Path:
        """Generate A3M file with pre-computed MSAs."""
        output_path = output_dir / f"{target_features.target_id}.a3m"

        # For complexes, we need to pair MSAs
        # For now, concatenate with gap characters

        with open(output_path, "w") as f:
            # Write query sequence first
            all_seqs = []
            for chain in target_features.chains:
                all_seqs.append(chain.sequence)

            combined_query = ":".join(all_seqs)
            f.write(f">{target_features.target_id}\n{combined_query}\n")

            # Write MSA sequences
            # Load each chain's MSA and combine
            chain_msas = []
            for chain in target_features.chains:
                if chain.msa_paths:
                    combined_chain_msas = []
                    seen = set()
                    for path in chain.msa_paths:
                        if path.exists():
                            for _, seq in self.load_msa(path):
                                s_clean = seq.replace("-", "").replace(".", "").upper()
                                if s_clean not in seen and s_clean != chain.sequence:
                                    seen.add(s_clean)
                                    combined_chain_msas.append((_, s_clean))
                    chain_msas.append(combined_chain_msas)
                else:
                    chain_msas.append([])

            # Find minimum MSA depth (for pairing)
            min_depth = min(len(m) for m in chain_msas) if chain_msas else 0

            for i in range(min(min_depth, 2048)):
                combined_seq = []
                for chain_idx, chain in enumerate(target_features.chains):
                    if i < len(chain_msas[chain_idx]):
                        _, seq = chain_msas[chain_idx][i]
                        # Pad or truncate to match chain length
                        if len(seq) < len(chain.sequence):
                            seq = seq + "-" * (len(chain.sequence) - len(seq))
                        else:
                            seq = seq[: len(chain.sequence)]
                        combined_seq.append(seq)
                    else:
                        # Gap for this chain
                        combined_seq.append("-" * len(chain.sequence))

                f.write(f">msa_{i}\n{':'.join(combined_seq)}\n")

        logger.info(f"Generated ColabFold A3M: {output_path}")
        return output_path


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
    parser = argparse.ArgumentParser(description="Generate ColabFold input")
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

    generator = ColabFoldInputGenerator(config)
    output_path = generator.generate(target_features, args.output_dir / args.target_id)

    if generator.validate_output(output_path):
        logger.info("Input generation successful")
        return 0
    else:
        logger.error("Input generation failed")
        return 1


if __name__ == "__main__":
    exit(main())
