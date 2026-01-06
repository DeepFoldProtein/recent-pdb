#!/usr/bin/env python3
"""
Target Registry Script
======================
Filters PDB entries based on benchmark criteria and generates sequence hash mappings.

Features:
- Resolution, date, and obsolete status filtering
- SHA256 sequence hashing for deduplication
- Chain-to-hash mapping for homo-oligomer optimization
- Output: target lists with sequence hash associations
"""

import argparse
import hashlib
import json
import logging
import sqlite3
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional

import yaml

try:
    import gemmi
except ImportError:
    gemmi = None

try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


@dataclass
class ChainInfo:
    """Information about a single chain in a PDB entry."""

    chain_id: str
    sequence: str
    seq_hash: str
    entity_id: Optional[str] = None


@dataclass
class TargetEntry:
    """A benchmark target with associated chain information."""

    pdb_id: str
    release_date: str
    resolution: Optional[float]
    chains: list[ChainInfo]
    file_path: str

    @property
    def unique_seq_hashes(self) -> set[str]:
        return {c.seq_hash for c in self.chains}


def compute_seq_hash(sequence: str) -> str:
    """Compute SHA256 hash of a protein sequence."""
    # Normalize: uppercase, remove whitespace
    normalized = sequence.upper().replace(" ", "").replace("\n", "")
    return hashlib.sha256(normalized.encode()).hexdigest()[:16]


def extract_sequences_gemmi(cif_path: Path) -> list[ChainInfo]:
    """Extract chain sequences from mmCIF using gemmi."""
    if gemmi is None:
        raise ImportError("gemmi is required for sequence extraction")

    chains = []

    try:
        structure = gemmi.read_structure(str(cif_path))

        for model in structure:
            for chain in model:
                # Get polymer sequence
                polymer = chain.get_polymer()
                if not polymer:
                    continue

                # Extract one-letter sequence
                sequence = polymer.make_one_letter_sequence()
                if len(sequence) < 10:  # Skip very short chains
                    continue

                seq_hash = compute_seq_hash(sequence)

                chains.append(
                    ChainInfo(
                        chain_id=chain.name,
                        sequence=sequence,
                        seq_hash=seq_hash,
                        # entity_id=chain.subchain, # invalid in gemmi 0.7.4
                    )
                )

            break  # Only process first model

    except Exception as e:
        logger.error(f"Failed to extract sequences from {cif_path}: {e}")

    return chains


def extract_sequences_db(seqres_str: str) -> list[ChainInfo]:
    """Extract chain sequences from stored SEQRES JSON."""
    import json

    data = json.loads(seqres_str)
    entities = data.get("entities", {})
    chains_map = data.get("chains", {})

    chains = []
    for asym_id, entity_id in chains_map.items():
        sequence = entities.get(entity_id)
        if sequence:
            chains.append(
                ChainInfo(
                    chain_id=asym_id,
                    sequence=sequence,
                    seq_hash=compute_seq_hash(sequence),
                    entity_id=entity_id,
                )
            )
    return chains


def filter_targets(
    conn: sqlite3.Connection,
    config: dict,
) -> list[TargetEntry]:
    """Filter PDB entries based on benchmark criteria."""
    filters = config.get("filters", {})

    query = """
        SELECT pdb_id, release_date, revision_date, resolution, file_path, seqres
        FROM pdb_entries
        WHERE 1=1
    """
    params = []

    # Resolution filter
    if "resolution_max" in filters:
        query += " AND (resolution IS NULL OR resolution <= ?)"
        params.append(filters["resolution_max"])

    # Date filter
    if "release_date_after" in filters:
        query += " AND release_date >= ?"
        params.append(filters["release_date_after"])

    # Obsolete filter
    if filters.get("exclude_obsolete", True):
        query += " AND is_obsolete = 0"

    query += " ORDER BY release_date DESC"

    cursor = conn.execute(query, params)

    targets = []
    min_len = filters.get("min_seq_length", 30)
    max_len = filters.get("max_seq_length", 2048)

    for row in cursor:
        pdb_id, release_date, revision_date, resolution, file_path, seqres = row

        if not file_path or not Path(file_path).exists():
            continue

        # Extract chain sequences
        if seqres:
            chains = extract_sequences_db(seqres)
        else:
            chains = extract_sequences_gemmi(Path(file_path))

        # Filter by sequence length
        chains = [c for c in chains if min_len <= len(c.sequence) <= max_len]

        if not chains:
            continue

        targets.append(
            TargetEntry(
                pdb_id=pdb_id,
                release_date=release_date,
                resolution=resolution,
                chains=chains,
                file_path=str(file_path),
            )
        )

    return targets


def load_custom_targets(
    conn: sqlite3.Connection,
    custom_file: Path,
    config: dict,
) -> list[TargetEntry]:
    """Load targets from a custom PDB ID file.

    Args:
        conn: SQLite database connection
        custom_file: Path to text file with one PDB ID per line
        config: Configuration dict (for sequence length filters)

    Returns:
        List of TargetEntry objects for the specified PDB IDs
    """
    filters = config.get("filters", {})
    min_len = filters.get("min_seq_length", 30)
    max_len = filters.get("max_seq_length", 2048)

    # Read PDB IDs from file (skip comments and empty lines)
    pdb_ids = []
    with open(custom_file) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                # Normalize: lowercase, remove extensions
                pdb_id = line.lower().split(".")[0]
                pdb_ids.append(pdb_id)

    logger.info(f"Loaded {len(pdb_ids)} PDB IDs from {custom_file}")

    # Query database for each PDB ID
    targets = []
    not_found = []

    for pdb_id in pdb_ids:
        cursor = conn.execute(
            """
            SELECT pdb_id, release_date, resolution, file_path, seqres
            FROM pdb_entries
            WHERE LOWER(pdb_id) = ?
            """,
            (pdb_id,),
        )
        row = cursor.fetchone()

        if not row:
            not_found.append(pdb_id)
            continue

        pdb_id_db, release_date, resolution, file_path, seqres = row

        if not file_path or not Path(file_path).exists():
            logger.warning(f"File not found for {pdb_id}: {file_path}")
            continue

        # Extract chain sequences
        if seqres:
            chains = extract_sequences_db(seqres)
        else:
            chains = extract_sequences_gemmi(Path(file_path))
        chains = [c for c in chains if min_len <= len(c.sequence) <= max_len]

        if not chains:
            logger.warning(f"No valid chains for {pdb_id}")
            continue

        targets.append(
            TargetEntry(
                pdb_id=pdb_id_db,
                release_date=release_date or "",
                resolution=resolution,
                chains=chains,
                file_path=str(file_path),
            )
        )

    if not_found:
        logger.warning(
            f"{len(not_found)} PDB IDs not found in registry: {not_found[:10]}..."
        )

    return targets


def build_hash_mapping(targets: list[TargetEntry]) -> dict[str, list[str]]:
    """Build mapping from sequence hash to target IDs."""
    mapping = {}

    for target in targets:
        for chain in target.chains:
            if chain.seq_hash not in mapping:
                mapping[chain.seq_hash] = []
            if target.pdb_id not in mapping[chain.seq_hash]:
                mapping[chain.seq_hash].append(target.pdb_id)

    return mapping


def write_outputs(
    targets: list[TargetEntry],
    hash_mapping: dict[str, list[str]],
    output_dir: Path,
):
    """Write target lists and mappings to output files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Write target IDs
    targets_file = output_dir / "targets.txt"
    with open(targets_file, "w") as f:
        for target in targets:
            f.write(f"{target.pdb_id}\n")
    logger.info(f"Wrote {len(targets)} targets to {targets_file}")

    # Write unique sequence hashes
    all_hashes = set()
    for target in targets:
        all_hashes.update(target.unique_seq_hashes)

    hashes_file = output_dir / "seq_hashes.txt"
    with open(hashes_file, "w") as f:
        for h in sorted(all_hashes):
            f.write(f"{h}\n")
    logger.info(f"Wrote {len(all_hashes)} unique sequence hashes to {hashes_file}")

    # Write hash-to-targets mapping (for MSA sharing)
    mapping_file = output_dir / "seq_to_targets.json"
    with open(mapping_file, "w") as f:
        json.dump(hash_mapping, f, indent=2)
    logger.info(f"Wrote hash mapping to {mapping_file}")

    # Write detailed target info (for debugging/analysis)
    details_file = output_dir / "target_details.json"
    details = []
    for target in targets:
        details.append(
            {
                "pdb_id": target.pdb_id,
                "release_date": target.release_date,
                "resolution": target.resolution,
                "num_chains": len(target.chains),
                "unique_sequences": len(target.unique_seq_hashes),
                "file_path": target.file_path,
                "chains": [asdict(c) for c in target.chains],
            }
        )
    with open(details_file, "w") as f:
        json.dump(details, f, indent=2)
    logger.info(f"Wrote target details to {details_file}")

    # Write ground truth mapping
    gt_file = output_dir / "ground_truths.json"
    gt_mapping = {t.pdb_id: t.file_path for t in targets}
    with open(gt_file, "w") as f:
        json.dump(gt_mapping, f, indent=2)
    logger.info(f"Wrote ground truth mapping to {gt_file}")

    # Write sequences for MSA search (one file per hash)
    seqs_dir = output_dir / "sequences"
    seqs_dir.mkdir(exist_ok=True)

    hash_to_seq = {}
    for target in targets:
        for chain in target.chains:
            if chain.seq_hash not in hash_to_seq:
                hash_to_seq[chain.seq_hash] = chain.sequence

    for seq_hash, sequence in hash_to_seq.items():
        seq_file = seqs_dir / f"{seq_hash}.fasta"
        with open(seq_file, "w") as f:
            f.write(f">{seq_hash}\n{sequence}\n")

    logger.info(f"Wrote {len(hash_to_seq)} sequence files to {seqs_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Filter PDB entries and generate target registry"
    )
    parser.add_argument(
        "--registry-db",
        type=Path,
        required=True,
        help="Path to SQLite registry database",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config.yaml"),
        help="Path to configuration file",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for target lists",
    )

    args = parser.parse_args()

    # Load config
    with open(args.config) as f:
        config = yaml.safe_load(f)

    # Connect to registry
    if not args.registry_db.exists():
        logger.error(f"Registry database not found: {args.registry_db}")
        logger.error("Run sync_pdb.py first to create the registry")
        return 1

    conn = sqlite3.connect(args.registry_db)

    try:
        # Check for custom targets file
        custom_file = config.get("filters", {}).get("custom_targets_file")

        if custom_file:
            custom_path = Path(custom_file)
            if not custom_path.exists():
                logger.error(f"Custom targets file not found: {custom_path}")
                return 1
            logger.info(f"Using custom targets from: {custom_path}")
            targets = load_custom_targets(conn, custom_path, config)
        else:
            # Filter targets using standard criteria
            logger.info("Filtering targets based on benchmark criteria...")
            targets = filter_targets(conn, config)

        logger.info(f"Found {len(targets)} qualifying targets")

        # Build hash mapping
        logger.info("Building sequence hash mappings...")
        hash_mapping = build_hash_mapping(targets)

        # Compute deduplication stats
        total_chains = sum(len(t.chains) for t in targets)
        unique_seqs = len(hash_mapping)
        logger.info(
            f"Deduplication: {total_chains} chains -> {unique_seqs} unique sequences "
            f"({100 * (1 - unique_seqs / max(total_chains, 1)):.1f}% reduction)"
        )

        # Write outputs
        write_outputs(targets, hash_mapping, args.output_dir)

    finally:
        conn.close()

    logger.info("Done!")
    return 0


if __name__ == "__main__":
    exit(main())
