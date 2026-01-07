#!/usr/bin/env python3
"""
MSA Search Controller
=====================
Runs MSA search for a single sequence hash using configured tools.

Features:
- Idempotent execution (skips if results exist)
- Specific tool and database targeting
- Result merging and deduplication
"""

import argparse
import logging
import subprocess
import tempfile
from pathlib import Path


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


class MSASearcher:
    """Base class for MSA search tools."""

    name: str = "base"

    def __init__(
        self,
        config: dict,
        container: Path = None,
        threads: int = 1,
        extra_binds: list[str] = None,
    ):
        self.config = config
        self.container = container
        self.threads = threads
        self.extra_binds = extra_binds or []

    def _get_apptainer_cmd(self, cmd: list, paths: list = None) -> list:
        """Wrap command with apptainer exec and necessary binds."""
        if not self.container:
            return cmd

        binds = set(self.extra_binds)
        if paths:
            for p in paths:
                if p:
                    p_path = Path(p)
                    if p_path.is_absolute():
                        # Bind the directory to ensure accessibility
                        # If it's a file, bind parent. If dir, bind itself or parent.
                        # Bind parent is generally safer for MMseqs databases (multi-file)
                        binds.add(str(p_path.parent))

        app_cmd = ["apptainer", "exec"]
        for b in sorted(binds):
            app_cmd.extend(["--bind", f"{b}:{b}"])

        app_cmd.append(str(self.container))
        return app_cmd + cmd

    def search(
        self,
        query_fasta: Path,
        output_path: Path,
        database: str,
    ) -> bool:
        """Run MSA search and save to output path."""
        raise NotImplementedError


class JackhmmerSearcher(MSASearcher):
    """JackHMMER MSA search."""

    name = "jackhmmer"

    def search(
        self,
        query_fasta: Path,
        output_path: Path,
        database: str,
    ) -> bool:
        params = self.config.get("msa_params", {}).get("jackhmmer", {})
        n_iter = params.get("n_iterations", 1)
        e_value = params.get("e_value", 1e-4)
        cpus = self.threads  # Use threads from init

        if output_path.exists():
            logger.info(f"Skipping jackhmmer (output exists): {output_path}")
            return True

        output_sto = output_path.with_suffix(".sto")

        cmd = [
            "jackhmmer",
            "-N",
            str(n_iter),
            "-E",
            str(e_value),
            "--cpu",
            str(cpus),
            "-A",
            str(output_sto),
            str(query_fasta),
            database,
        ]

        logger.info(f"Running: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600 * 4,  # 4 hour timeout
            )

            if result.returncode != 0:
                logger.error(f"jackhmmer failed: {result.stderr}")
                return False

            # Convert STO to A3M
            self._sto_to_a3m(output_sto, output_path)
            # Cleanup STO
            if output_sto.exists():
                output_sto.unlink()
            return True

        except subprocess.TimeoutExpired:
            logger.error("jackhmmer timed out")
            return False
        except FileNotFoundError:
            logger.error("jackhmmer not found in PATH")
            return False

    def _sto_to_a3m(self, sto_path: Path, a3m_path: Path):
        """Convert Stockholm to A3M format."""
        # Simple conversion - in production use proper tool
        sequences = {}
        with open(sto_path) as f:
            for line in f:
                if line.startswith("#") or line.startswith("//") or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    name, seq = parts[0], parts[1]
                    if name not in sequences:
                        sequences[name] = ""
                    sequences[name] += seq

        with open(a3m_path, "w") as f:
            for name, seq in sequences.items():
                # Remove gaps for A3M format
                clean_seq = seq.replace(".", "").replace("-", "")
                if clean_seq:
                    f.write(f">{name}\n{seq}\n")


class HHblitsSearcher(MSASearcher):
    """HHblits MSA search."""

    name = "hhblits"

    def search(
        self,
        query_fasta: Path,
        output_path: Path,
        database: str,
    ) -> bool:
        # Note: config is now passed from msa.tools.{tool_name}.params usually
        # But we kept partial compatibility or rely on full structure
        # Let's check config structure passed to init.
        # It's the full config object.
        params = (
            self.config.get("msa", {})
            .get("tools", {})
            .get("hhblits", {})
            .get("params", {})
        )

        # Fallback if config structure varies
        if not params:
            params = self.config.get("msa_params", {}).get("hhblits", {})

        n_iter = params.get("n_iterations", 3)
        n_iter = params.get("n_iterations", 3)
        # e_pre removed (unused)
        e_value = params.get("e_value", 0.001)
        cpus = self.threads  # Use threads from init

        if output_path.exists():
            logger.info(f"Skipping hhblits (output exists): {output_path}")
            return True

        cmd = [
            "hhblits",
            "-i",
            str(query_fasta),
            "-o",
            "/dev/null",
            "-oa3m",
            str(output_path),
            "-d",
            database,
            "-n",
            str(n_iter),
            "-e",
            str(e_value),
            "-cpu",
            str(cpus),
        ]

        if self.container:
            cmd = self._get_apptainer_cmd(
                cmd, paths=[query_fasta, output_path, database]
            )

        logger.info(f"Running: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600 * 4,
            )

            if result.returncode != 0:
                logger.error(f"hhblits failed: {result.stderr}")
                return False

            return True

        except subprocess.TimeoutExpired:
            logger.error("hhblits timed out")
            return False
        except FileNotFoundError:
            logger.error("hhblits not found in PATH")
            return False


class MMseqs2Searcher(MSASearcher):
    """MMseqs2 MSA search."""

    name = "mmseqs2"

    def search(
        self,
        query_fasta: Path,
        output_path: Path,
        database: str,
    ) -> bool:
        # Config params
        params = (
            self.config.get("msa", {})
            .get("tools", {})
            .get("mmseqs2", {})
            .get("params", {})
        )
        params = params if params else {}

        if output_path.exists():
            logger.info(f"Skipping mmseqs2 (output exists): {output_path}")
            return True

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            query_db = tmp / "query_db"
            result_db = tmp / "result_db"

            # 1. Create query database
            mmseqs_bin = "/usr/local/bin/entrypoint" if self.container else "mmseqs"
            createdb_cmd = [
                mmseqs_bin,
                "createdb",
                str(query_fasta),
                str(query_db),
            ]
            if self.container:
                createdb_cmd = self._get_apptainer_cmd(
                    createdb_cmd, paths=[query_fasta, query_db]
                )

            logger.info(f"Running createdb: {' '.join(createdb_cmd)}")
            try:
                subprocess.run(createdb_cmd, capture_output=True, check=True, text=True)
            except subprocess.CalledProcessError as e:
                logger.error(f"mmseqs createdb failed: {e.stderr}")
                return False
            except FileNotFoundError:
                logger.error("mmseqs not found in PATH")
                return False

            # 2. MMseqs search
            search_cmd = [
                mmseqs_bin,
                "search",
                str(query_db),
                database,
                str(result_db),
                str(tmp),
                "--num-iterations",
                "3",  # Usually 3 is good for sensitivity
                "--threads",
                str(self.threads),
                *(item for k, v in params.items() for item in (f"--{k}", str(v))),
            ]

            if self.container:
                search_cmd = self._get_apptainer_cmd(
                    search_cmd, paths=[query_db, database, tmp]
                )

            logger.info(f"Running search: {' '.join(search_cmd)}")

            try:
                result = subprocess.run(
                    search_cmd,
                    capture_output=True,
                    text=True,
                    timeout=3600 * 2,
                )

                if result.returncode != 0:
                    logger.error(f"mmseqs search failed: {result.stderr}")
                    return False

                # 3. Result to MSA (A3M)
                msa_cmd = [
                    mmseqs_bin,
                    "result2msa",
                    str(query_db),
                    database,
                    str(result_db),
                    str(output_path),
                    "--threads",
                    str(self.threads),
                ]

                if self.container:
                    msa_cmd = self._get_apptainer_cmd(
                        msa_cmd, paths=[query_db, database, output_path]
                    )

                logger.info(f"Running result2msa: {' '.join(msa_cmd)}")

                result = subprocess.run(
                    msa_cmd,
                    capture_output=True,
                    text=True,
                    timeout=600,
                )

                if result.returncode != 0:
                    logger.error(f"mmseqs result2msa failed: {result.stderr}")
                    return False

                return True

            except subprocess.TimeoutExpired:
                logger.error("mmseqs2 timed out")
                return False
            except FileNotFoundError:
                logger.error("mmseqs not found in PATH")
                return False


def run_msa_search(
    seq_hash: str,
    tool_name: str,
    db_name: str,
    db_version: str,
    output_path: Path,
    cache_dir: Path,
    config: dict,
    container: Path = None,
    threads: int = 1,
    extra_binds: list[str] = None,
):
    """Run MSA search for a single sequence hash, tool, and database."""
    # Setup output directory
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Find query sequence file
    # Legacy cache dir for fasta might be here?
    # Actually Snakefile provides query input, but here we scan if not passed?
    # Snakefile only passed seq-hash. It presumed query file location.
    # We need to find the query file.

    query_fasta = None
    target_lists = Path(config["paths"]["target_lists"])
    potential_locations = [
        target_lists / "sequences" / f"{seq_hash}.fasta",  # Standard location
        cache_dir / seq_hash / "query.fasta",  # Legacy location
    ]

    for loc in potential_locations:
        if loc.exists():
            query_fasta = loc
            break

    if query_fasta is None:
        logger.error(f"Query sequence not found for hash: {seq_hash}")
        return False

    # Initialize searchers
    searchers = {
        "jackhmmer": JackhmmerSearcher(config, container, threads),
        "hhblits": HHblitsSearcher(config, container, threads),
        "mmseqs2": MMseqs2Searcher(config, container, threads),
    }

    if tool_name not in searchers:
        logger.error(f"Unknown tool: {tool_name}")
        return False

    searcher = searchers[tool_name]

    # Find database path from config
    db_path = None
    tools_config = config.get("msa", {}).get("tools", {})
    tool_config = tools_config.get(tool_name, {})
    databases = tool_config.get("databases", [])

    found = False
    for db in databases:
        if db["name"] == db_name and str(db.get("version", "")) == str(db_version):
            db_path = db["path"]
            found = True
            break

    if not found:
        # Fallback: check if db_name matches one entry without version check if version is empty
        # But for correctness let's be strict or lenient based on request
        # If db_version is "None" or empty string, maybe match matches with empty version?
        # Let's try to match by name only if exact match failed and db_version is empy
        if not db_version:
            for db in databases:
                if db["name"] == db_name:
                    db_path = db["path"]
                    found = True
                    break

    if not db_path:
        logger.error(
            f"Database not found in config: tool={tool_name}, db={db_name}, version={db_version}"
        )
        return False

    logger.info(
        f"Starting {tool_name} search on {db_name} ({db_version}) for {seq_hash}"
    )
    success = searcher.search(query_fasta, output_path, db_path)

    if success:
        logger.info(f"MSA saved to: {output_path}")
    else:
        logger.error("MSA search failed")

    return success


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
    parser = argparse.ArgumentParser(description="Run MSA search for a sequence hash")
    parser.add_argument("--seq-hash", type=str, required=True, help="Sequence hash")
    parser.add_argument("--tool", type=str, required=True, help="MSA tool name")
    parser.add_argument("--db-name", type=str, required=True, help="Database name")
    parser.add_argument("--db-version", type=str, default="", help="Database version")
    parser.add_argument(
        "--output", type=Path, required=True, help="Output A3M file path"
    )
    parser.add_argument(
        "--cache-dir", type=Path, required=True, help="Cache directory root"
    )
    parser.add_argument("--container", type=Path, help="Path to Apptainer container")
    parser.add_argument("--threads", type=int, default=1, help="Number of CPU threads")
    parser.add_argument(
        "--extra-binds",
        type=str,
        help="Comma-separated list of extra bind paths for Apptainer",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config.yaml"),
        help="Path to configuration file",
    )

    args = parser.parse_args()

    # Load config
    config = _load_config(args.config)

    # Parse extra binds
    extra_binds = []
    if args.extra_binds:
        extra_binds = [b.strip() for b in args.extra_binds.split(",")]

    # Run search
    success = run_msa_search(
        args.seq_hash,
        args.tool,
        args.db_name,
        args.db_version,
        args.output,
        args.cache_dir,
        config,
        args.container,
        args.threads,
        extra_binds,
    )

    return 0 if success else 1


if __name__ == "__main__":
    exit(main())
