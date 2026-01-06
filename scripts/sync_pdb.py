#!/usr/bin/env python3
"""
PDB Sync & Audit Script
=======================
Synchronizes local PDB mirror with remote source and updates metadata registry.

Features:
- Incremental sync using rsync (only changed files)
- mmCIF header parsing for metadata extraction
- SQLite registry with revision tracking
- Obsolete entry detection
"""

import argparse
import logging
import sqlite3
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterator, Optional

import yaml

try:
    import gemmi
except ImportError:
    gemmi = None

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


@dataclass
class PDBMetadata:
    """Metadata extracted from a PDB mmCIF file."""

    pdb_id: str
    release_date: Optional[str]
    revision_date: Optional[str]
    resolution: Optional[float]
    method: Optional[str]
    title: Optional[str]
    is_obsolete: bool
    file_path: str
    file_mtime: float
    seqres: Optional[str] = None


def init_database(db_path: Path) -> sqlite3.Connection:
    """Initialize SQLite database with schema."""
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS pdb_entries (
            pdb_id TEXT PRIMARY KEY,
            release_date TEXT,
            revision_date TEXT,
            resolution REAL,
            method TEXT,
            title TEXT,
            is_obsolete INTEGER DEFAULT 0,
            file_path TEXT,
            file_mtime REAL,
            seqres TEXT,
            last_updated TEXT
        )
    """)
    conn.execute("""
        CREATE INDEX IF NOT EXISTS idx_release_date ON pdb_entries(release_date)
    """)
    conn.execute("""
        CREATE INDEX IF NOT EXISTS idx_resolution ON pdb_entries(resolution)
    """)
    conn.commit()
    return conn


def parse_mmcif_metadata(cif_path: Path) -> Optional[PDBMetadata]:
    """Extract metadata from mmCIF file using gemmi."""
    if gemmi is None:
        logger.warning("gemmi not installed, skipping metadata extraction")
        return None

    try:
        doc = gemmi.cif.read(str(cif_path))
        block = doc.sole_block()

        pdb_id = block.find_value("_entry.id")
        if pdb_id:
            pdb_id = pdb_id.strip().upper()

        # Release date
        release_date = None
        for tag in [
            "_pdbx_database_status.recvd_initial_deposition_date",
            "_database_PDB_rev.date_original",
        ]:
            val = block.find_value(tag)
            if val and val != "?":
                release_date = val.strip()
                break

        # Revision date
        revision_date = None
        rev_table = block.find(["_pdbx_audit_revision_history.revision_date"])
        if rev_table:
            dates = [row[0] for row in rev_table if row[0] != "?"]
            if dates:
                revision_date = max(dates)

        # Resolution
        resolution = None
        for tag in ["_refine.ls_d_res_high", "_em_3d_reconstruction.resolution"]:
            val = block.find_value(tag)
            if val and val not in ("?", "."):
                try:
                    resolution = float(val)
                    break
                except ValueError:
                    continue

        # Experimental method
        method = block.find_value("_exptl.method")
        if method:
            method = method.strip().strip("'\"")

        # Title
        title = block.find_value("_struct.title")
        if title:
            title = title.strip().strip("'\"")[:500]  # Truncate long titles

        # Obsolete check
        is_obsolete = False
        status = block.find_value("_pdbx_database_status.status_code")
        if status and "OBS" in status.upper():
            is_obsolete = True

        # SEQRES extraction
        seqres_dict = {}
        poly_seq = block.get_mmcif_category("_entity_poly_seq.")
        if "entity_id" in poly_seq and "mon_id" in poly_seq:
            from collections import defaultdict

            eid_to_res = defaultdict(list)
            for eid, mon in zip(poly_seq["entity_id"], poly_seq["mon_id"]):
                one_letter = gemmi.one_letter_code([mon])
                if not one_letter or one_letter == " ":
                    one_letter = "X"
                eid_to_res[eid].append(one_letter)
            for eid, res_list in eid_to_res.items():
                seqres_dict[eid] = "".join(res_list)

        # Mapping asym_id to entity_id
        asym_to_entity = {}
        asym = block.get_mmcif_category("_struct_asym.")
        if "id" in asym and "entity_id" in asym:
            for aid, eid in zip(asym["id"], asym["entity_id"]):
                if eid in seqres_dict:
                    asym_to_entity[aid] = eid

        import json

        seqres_json = json.dumps({"entities": seqres_dict, "chains": asym_to_entity})

        return PDBMetadata(
            pdb_id=pdb_id or cif_path.stem.upper(),
            release_date=release_date,
            revision_date=revision_date,
            resolution=resolution,
            method=method,
            title=title,
            is_obsolete=is_obsolete,
            file_path=str(cif_path),
            file_mtime=cif_path.stat().st_mtime,
            seqres=seqres_json,
        )
    except Exception as e:
        logger.error(f"Failed to parse {cif_path}: {e}")
        return None


def scan_pdb_directory(pdb_dir: Path) -> Iterator[Path]:
    """Recursively find all mmCIF files in PDB directory."""
    for pattern in ["**/*.cif", "**/*.cif.gz"]:
        yield from pdb_dir.glob(pattern)


def sync_pdb_mirror(source: str, dest: Path, dry_run: bool = False) -> bool:
    """Run rsync to synchronize PDB mirror."""
    cmd = [
        "rsync",
        "-avz",
        "--delete",
        "--include=*/",
        "--include=*.cif.gz",
        "--exclude=*",
    ]

    if dry_run:
        cmd.append("--dry-run")

    cmd.extend([source, str(dest)])

    logger.info(f"Running rsync: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logger.error(f"rsync failed: {result.stderr}")
        return False

    logger.info("rsync completed successfully")
    return True


def update_registry(
    conn: sqlite3.Connection, pdb_dir: Path, force_rescan: bool = False
) -> tuple[int, int, int]:
    """Update registry with new/modified PDB entries."""
    added = 0
    updated = 0
    skipped = 0

    cursor = conn.cursor()

    for cif_path in scan_pdb_directory(pdb_dir):
        pdb_id = cif_path.stem.split(".")[0].upper()
        file_mtime = cif_path.stat().st_mtime

        row = None
        # Check if entry exists and is up-to-date
        if not force_rescan:
            cursor.execute(
                "SELECT file_mtime FROM pdb_entries WHERE pdb_id = ?", (pdb_id,)
            )
            row = cursor.fetchone()
            if row and abs(row[0] - file_mtime) < 1.0:
                skipped += 1
                continue

        # Parse metadata
        metadata = parse_mmcif_metadata(cif_path)
        if metadata is None:
            continue

        # Upsert entry
        cursor.execute(
            """
            INSERT OR REPLACE INTO pdb_entries
            (pdb_id, release_date, revision_date, resolution, method, title,
             is_obsolete, file_path, file_mtime, seqres, last_updated)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
            (
                metadata.pdb_id,
                metadata.release_date,
                metadata.revision_date,
                metadata.resolution,
                metadata.method,
                metadata.title,
                1 if metadata.is_obsolete else 0,
                metadata.file_path,
                metadata.file_mtime,
                metadata.seqres,
                datetime.now().isoformat(),
            ),
        )

        if row:
            updated += 1
        else:
            added += 1

    conn.commit()
    return added, updated, skipped


def main():
    parser = argparse.ArgumentParser(
        description="Synchronize PDB mirror and update metadata registry"
    )
    parser.add_argument(
        "--pdb-master", type=Path, help="Path to central PDB mirror directory"
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
        "--rsync-source",
        type=str,
        default="rsync.rcsb.org::ftp_data/structures/divided/mmCIF/",
        help="rsync source URL for PDB mirror",
    )
    parser.add_argument(
        "--skip-sync", action="store_true", help="Skip rsync and only update registry"
    )
    parser.add_argument(
        "--force-rescan", action="store_true", help="Force rescan of all PDB files"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without making changes",
    )

    args = parser.parse_args()

    # Load config if pdb_master not specified
    if args.pdb_master is None:
        with open(args.config) as f:
            config = yaml.safe_load(f)
        args.pdb_master = Path(config["paths"]["pdb_master"])

    # Ensure directories exist
    args.registry_db.parent.mkdir(parents=True, exist_ok=True)

    # Run rsync if not skipped
    if not args.skip_sync:
        logger.info("Starting PDB mirror synchronization...")
        if not sync_pdb_mirror(args.rsync_source, args.pdb_master, args.dry_run):
            sys.exit(1)

    # Update registry
    logger.info("Updating metadata registry...")
    conn = init_database(args.registry_db)

    try:
        added, updated, skipped = update_registry(
            conn, args.pdb_master, args.force_rescan
        )
        logger.info(
            f"Registry update complete: {added} added, {updated} updated, "
            f"{skipped} unchanged"
        )
    finally:
        conn.close()

    logger.info("Done!")


if __name__ == "__main__":
    main()
