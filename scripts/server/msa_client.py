#!/usr/bin/env python3
"""
MSA Server Client
=================
Client for the in-house MSA server. Compatible with ColabFold API protocol.
Can be used to submit MSA jobs and retrieve results.
"""

import argparse
import logging
import time
from pathlib import Path
from typing import Optional

try:
    import requests
except ImportError:
    print("requests not installed. Run: pip install requests")
    raise

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


class MSAClient:
    """Client for MSA server."""

    def __init__(
        self,
        server_url: str = "http://localhost:8000",
        timeout: int = 3600,
        poll_interval: int = 5,
    ):
        self.server_url = server_url.rstrip("/")
        self.timeout = timeout
        self.poll_interval = poll_interval

    def search(
        self,
        sequence: str,
        databases: list[str] = None,
        mode: str = "all",
    ) -> Optional[str]:
        """
        Submit MSA search and wait for results.

        Args:
            sequence: Protein sequence (FASTA format or raw)
            databases: List of database names to search
            mode: Search mode ("all", "unpaired", "paired")

        Returns:
            A3M format MSA string, or None on failure
        """
        if databases is None:
            databases = ["uniref30", "colabfold_envdb"]

        # Submit job
        response = requests.post(
            f"{self.server_url}/ticket/msa",
            json={
                "query": sequence,
                "databases": databases,
                "mode": mode,
            },
            timeout=30,
        )
        response.raise_for_status()

        ticket = response.json()
        job_id = ticket["id"]
        logger.info(f"Submitted MSA job: {job_id}")

        # Poll for results
        start_time = time.time()

        while time.time() - start_time < self.timeout:
            response = requests.get(
                f"{self.server_url}/ticket/msa/{job_id}",
                timeout=30,
            )
            response.raise_for_status()

            status = response.json()

            if status["status"] == "complete":
                logger.info(f"Job {job_id} complete")
                return status.get("a3m")

            elif status["status"] == "error":
                logger.error(f"Job {job_id} failed: {status.get('error')}")
                return None

            else:
                logger.debug(f"Job {job_id} status: {status['status']}")
                time.sleep(self.poll_interval)

        logger.error(f"Job {job_id} timed out")
        return None

    def search_sync(
        self,
        sequence: str,
        databases: list[str] = None,
    ) -> Optional[str]:
        """
        Direct synchronous MSA search.

        Args:
            sequence: Protein sequence
            databases: List of database names

        Returns:
            A3M format MSA string
        """
        if databases is None:
            databases = ["uniref30", "colabfold_envdb"]

        response = requests.post(
            f"{self.server_url}/result/msa",
            json={
                "query": sequence,
                "databases": databases,
            },
            timeout=self.timeout,
        )
        response.raise_for_status()

        return response.text

    def health_check(self) -> dict:
        """Check server health."""
        response = requests.get(
            f"{self.server_url}/health",
            timeout=10,
        )
        response.raise_for_status()
        return response.json()


def main():
    parser = argparse.ArgumentParser(description="MSA Server Client")
    parser.add_argument(
        "--server", default="http://localhost:8000", help="MSA server URL"
    )
    parser.add_argument("--input", type=Path, required=True, help="Input FASTA file")
    parser.add_argument("--output", type=Path, required=True, help="Output A3M file")
    parser.add_argument(
        "--databases",
        nargs="+",
        default=["uniref30", "colabfold_envdb"],
        help="Databases to search",
    )
    parser.add_argument(
        "--sync", action="store_true", help="Use synchronous (direct) search"
    )
    parser.add_argument("--timeout", type=int, default=3600, help="Timeout in seconds")

    args = parser.parse_args()

    # Read input sequence
    with open(args.input) as f:
        sequence = f.read()

    # Create client
    client = MSAClient(
        server_url=args.server,
        timeout=args.timeout,
    )

    # Check server health
    try:
        health = client.health_check()
        logger.info(f"Server status: {health}")
    except Exception as e:
        logger.error(f"Server not available: {e}")
        return 1

    # Run search
    if args.sync:
        a3m = client.search_sync(sequence, args.databases)
    else:
        a3m = client.search(sequence, args.databases)

    if a3m:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as f:
            f.write(a3m)
        logger.info(f"Saved MSA to: {args.output}")
        return 0
    else:
        logger.error("MSA search failed")
        return 1


if __name__ == "__main__":
    exit(main())
