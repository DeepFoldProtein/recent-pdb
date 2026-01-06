#!/usr/bin/env python3
"""
Test Setup Script
=================
Creates a minimal test dataset for pipeline verification.
Downloads a few recent PDB structures (gzipped) and prepares them for testing.
"""

import gzip
import hashlib
import json
import logging
import subprocess
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

# Test targets - small, recent proteins
TEST_TARGETS = [
    # PDB ID, Chain(s), Description
    ("8WGV", ["A"], "Small protein (~100 residues)"),
    ("8X4R", ["A"], "Single chain protein"),
    ("8Y0E", ["A", "B"], "Protein dimer"),
]


def compute_seq_hash(sequence: str) -> str:
    """Compute SHA256 hash of sequence."""
    normalized = sequence.upper().replace(" ", "").replace("\n", "")
    return hashlib.sha256(normalized.encode()).hexdigest()[:16]


def download_pdb_gz(pdb_id: str, output_dir: Path) -> Path:
    """Download PDB structure in gzipped mmCIF format."""
    output_file = output_dir / f"{pdb_id}.cif.gz"

    if output_file.exists():
        logger.info(f"  Already exists: {output_file}")
        return output_file

    url = f"https://files.rcsb.org/download/{pdb_id}.cif.gz"
    logger.info(f"  Downloading {url}...")

    try:
        subprocess.run(
            ["curl", "-sL", "-o", str(output_file), url],
            check=True,
            timeout=60,
        )
        return output_file
    except Exception as e:
        logger.error(f"  Failed to download {pdb_id}: {e}")
        return None


def read_cif_gz(cif_path: Path) -> str:
    """Read gzipped CIF file content."""
    with gzip.open(cif_path, "rt") as f:
        return f.read()


def extract_sequence_gemmi(cif_path: Path, chain_ids: list[str]) -> dict[str, str]:
    """Extract sequences from gzipped mmCIF using gemmi."""
    try:
        import gemmi
    except ImportError:
        logger.warning("gemmi not installed, using regex fallback")
        return extract_sequence_fallback(cif_path, chain_ids)

    sequences = {}

    try:
        # gemmi can read .gz files directly
        structure = gemmi.read_structure(str(cif_path))

        for model in structure:
            for chain in model:
                if chain.name in chain_ids or not chain_ids:
                    polymer = chain.get_polymer()
                    if polymer:
                        seq = polymer.make_one_letter_sequence()
                        if len(seq) >= 20:
                            sequences[chain.name] = seq
            break  # First model only
    except Exception as e:
        logger.error(f"  Failed to parse {cif_path}: {e}")

    return sequences


def extract_sequence_fallback(cif_path: Path, chain_ids: list[str]) -> dict[str, str]:
    """Extract sequences from gzipped mmCIF using regex (fallback)."""
    import re

    sequences = {}
    content = read_cif_gz(cif_path)

    # Look for _entity_poly.pdbx_seq_one_letter_code_can (canonical sequence)
    # This is more reliable than pdbx_seq_one_letter_code
    patterns = [
        r"_entity_poly\.pdbx_seq_one_letter_code_can\s*\n?['\"]?([A-Z\n\s]+)['\"]?",
        r"_entity_poly\.pdbx_seq_one_letter_code\s*\n?['\"]?([A-Z\n\s]+)['\"]?",
    ]

    for pattern in patterns:
        matches = re.findall(pattern, content, re.MULTILINE | re.DOTALL)
        if matches:
            for i, match in enumerate(matches):
                seq = match.replace("\n", "").replace(" ", "").strip()
                # Remove trailing semicolons or quotes
                seq = seq.rstrip(";'\"")
                if len(seq) >= 20:
                    # Use provided chain_ids or default to A, B, C...
                    if i < len(chain_ids):
                        chain_id = chain_ids[i]
                    else:
                        chain_id = chr(ord("A") + i)
                    sequences[chain_id] = seq
            break

    return sequences


def setup_test_data():
    """Set up test data directories and files."""
    base_dir = Path(".")
    target_lists = base_dir / "data" / "target_lists"
    ground_truths = base_dir / "data" / "ground_truths"
    sequences_dir = target_lists / "sequences"

    # Create directories
    target_lists.mkdir(parents=True, exist_ok=True)
    ground_truths.mkdir(parents=True, exist_ok=True)
    sequences_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Setting up test dataset")
    logger.info("=" * 60)

    targets = []
    seq_to_targets = {}
    all_hashes = set()
    target_details = []

    for pdb_id, chain_ids, description in TEST_TARGETS:
        logger.info(f"\nProcessing {pdb_id}: {description}")

        # Download gzipped structure
        cif_path = download_pdb_gz(pdb_id, ground_truths)
        if cif_path is None or not cif_path.exists():
            logger.warning(f"  Skipping {pdb_id}")
            continue

        # Extract sequences
        sequences = extract_sequence_gemmi(cif_path, chain_ids)
        if not sequences:
            logger.warning(f"  No sequences found for {pdb_id}")
            continue

        logger.info(f"  Found {len(sequences)} chain(s)")

        # Build target info
        chains = []
        for chain_id, seq in sequences.items():
            seq_hash = compute_seq_hash(seq)
            all_hashes.add(seq_hash)

            chains.append(
                {
                    "chain_id": chain_id,
                    "sequence": seq,
                    "seq_hash": seq_hash,
                }
            )

            # Write sequence FASTA
            seq_file = sequences_dir / f"{seq_hash}.fasta"
            with open(seq_file, "w") as f:
                f.write(f">{seq_hash}\n{seq}\n")

            # Update mapping
            if seq_hash not in seq_to_targets:
                seq_to_targets[seq_hash] = []
            if pdb_id not in seq_to_targets[seq_hash]:
                seq_to_targets[seq_hash].append(pdb_id)

            logger.info(
                f"    Chain {chain_id}: {len(seq)} residues, hash={seq_hash[:8]}..."
            )

        targets.append(pdb_id)
        target_details.append(
            {
                "pdb_id": pdb_id,
                "release_date": "2024-01-01",
                "resolution": 2.0,
                "num_chains": len(chains),
                "chains": chains,
            }
        )

    # Write output files
    logger.info("\n" + "=" * 60)
    logger.info("Writing output files")
    logger.info("=" * 60)

    # targets.txt
    targets_file = target_lists / "targets.txt"
    with open(targets_file, "w") as f:
        for t in targets:
            f.write(f"{t}\n")
    logger.info(f"  {targets_file}: {len(targets)} targets")

    # seq_hashes.txt
    hashes_file = target_lists / "seq_hashes.txt"
    with open(hashes_file, "w") as f:
        for h in sorted(all_hashes):
            f.write(f"{h}\n")
    logger.info(f"  {hashes_file}: {len(all_hashes)} unique sequences")

    # seq_to_targets.json
    mapping_file = target_lists / "seq_to_targets.json"
    with open(mapping_file, "w") as f:
        json.dump(seq_to_targets, f, indent=2)
    logger.info(f"  {mapping_file}")

    # target_details.json
    details_file = target_lists / "target_details.json"
    with open(details_file, "w") as f:
        json.dump(target_details, f, indent=2)
    logger.info(f"  {details_file}")

    logger.info("\n" + "=" * 60)
    logger.info("Test setup complete!")
    logger.info("=" * 60)
    logger.info(f"\nTargets: {targets}")
    logger.info(f"Unique sequences: {len(all_hashes)}")
    logger.info("\nNext steps:")
    logger.info("  1. Run: snakemake -n  (dry-run)")
    logger.info("  2. Run: snakemake run_all_msa --cores 4")

    return len(targets)


if __name__ == "__main__":
    n = setup_test_data()
    sys.exit(0 if n > 0 else 1)
