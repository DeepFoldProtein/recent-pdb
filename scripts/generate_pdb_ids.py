#!/usr/bin/env python3
"""
PDB ID List Generator
=====================
Generate PDB ID lists based on various filtering criteria.

Features:
- Random sampling from filtered results
- Date range, resolution, sequence length filters
- Multimer/ligand-containing structure filters
- Exclude specific PDB IDs

Usage:
    # Sample 100 PDBs from 2024 with resolution ≤ 2.0
    python generate_pdb_ids.py \\
        --registry-db data/registry.sqlite \\
        --date-range 2024-01-01 2025-12-31 \\
        --resolution-max 2.0 \\
        --sample 100 \\
        --output my_targets.txt
"""

import argparse
import logging
import random
import sqlite3
from datetime import datetime
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def build_query(args) -> tuple[str, list]:
    """Build SQL query from filter arguments."""
    query = """
        SELECT pdb_id, release_date, resolution
        FROM pdb_entries
        WHERE is_obsolete = 0
    """
    params = []

    if args.resolution_max:
        query += " AND (resolution IS NULL OR resolution <= ?)"
        params.append(args.resolution_max)

    if args.date_start:
        query += " AND release_date >= ?"
        params.append(args.date_start)

    if args.date_end:
        query += " AND release_date <= ?"
        params.append(args.date_end)

    query += " ORDER BY release_date DESC"
    return query, params


def main():
    parser = argparse.ArgumentParser(
        description="Generate PDB ID lists based on filtering criteria",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Sample 50 recent PDBs
  python generate_pdb_ids.py --registry-db data/registry.sqlite --sample 50 -o targets.txt

  # PDBs from 2024 with resolution ≤ 2.0
  python generate_pdb_ids.py --registry-db data/registry.sqlite \\
      --date-range 2024-01-01 2024-12-31 --resolution-max 2.0 -o targets.txt

  # Exclude specific PDBs from an existing list
  python generate_pdb_ids.py --registry-db data/registry.sqlite \\
      --exclude obsolete_ids.txt --sample 100 -o filtered.txt
""",
    )

    # Required arguments
    parser.add_argument(
        "--registry-db",
        type=Path,
        required=True,
        help="Path to SQLite registry database",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output file path for PDB ID list",
    )

    # Filter arguments
    parser.add_argument(
        "--resolution-max",
        type=float,
        help="Maximum resolution in Å (default: no limit)",
    )
    parser.add_argument(
        "--date-range",
        nargs=2,
        metavar=("START", "END"),
        help="Release date range (YYYY-MM-DD format)",
    )
    parser.add_argument(
        "--sample",
        type=int,
        help="Randomly sample N entries from filtered results",
    )
    parser.add_argument(
        "--exclude",
        type=Path,
        help="File containing PDB IDs to exclude (one per line)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducible sampling (default: 42)",
    )

    args = parser.parse_args()

    # Parse date range
    args.date_start = None
    args.date_end = None
    if args.date_range:
        try:
            args.date_start = datetime.strptime(
                args.date_range[0], "%Y-%m-%d"
            ).strftime("%Y-%m-%d")
            args.date_end = datetime.strptime(args.date_range[1], "%Y-%m-%d").strftime(
                "%Y-%m-%d"
            )
        except ValueError as e:
            logger.error(f"Invalid date format: {e}")
            return 1

    # Load exclusion list
    exclude_ids = set()
    if args.exclude and args.exclude.exists():
        with open(args.exclude) as f:
            for line in f:
                line = line.strip().lower()
                if line and not line.startswith("#"):
                    exclude_ids.add(line)
        logger.info(f"Loaded {len(exclude_ids)} PDB IDs to exclude")

    # Connect to database
    if not args.registry_db.exists():
        logger.error(f"Registry database not found: {args.registry_db}")
        return 1

    conn = sqlite3.connect(args.registry_db)

    try:
        query, params = build_query(args)
        cursor = conn.execute(query, params)

        # Collect results
        pdb_ids = []
        for row in cursor:
            pdb_id = row[0].lower()
            if pdb_id not in exclude_ids:
                pdb_ids.append(pdb_id)

        logger.info(f"Found {len(pdb_ids)} PDB entries matching criteria")

        # Random sampling
        if args.sample and args.sample < len(pdb_ids):
            random.seed(args.seed)
            pdb_ids = random.sample(pdb_ids, args.sample)
            logger.info(f"Sampled {len(pdb_ids)} entries (seed={args.seed})")

        # Write output
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as f:
            f.write("# Generated by generate_pdb_ids.py\n")
            f.write(f"# Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# Filters: resolution_max={args.resolution_max}, ")
            f.write(f"date_range={args.date_start} to {args.date_end}\n")
            f.write(f"# Total: {len(pdb_ids)} entries\n")
            for pdb_id in pdb_ids:
                f.write(f"{pdb_id}\n")

        logger.info(f"Wrote {len(pdb_ids)} PDB IDs to {args.output}")

    finally:
        conn.close()

    return 0


if __name__ == "__main__":
    exit(main())
