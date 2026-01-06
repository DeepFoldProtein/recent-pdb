#!/usr/bin/env python3
"""
MSA Merge Script
================
Merges multiple A3M files from different MSA tools.
Removes duplicate sequences while preserving order.
"""

import argparse
import logging
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


def parse_a3m(path: Path) -> list[tuple[str, str]]:
    """Parse A3M file into list of (header, sequence) tuples."""
    sequences = []

    if not path.exists():
        logger.warning(f"File not found: {path}")
        return sequences

    with open(path) as f:
        current_header = None
        current_seq = []

        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_header:
                    sequences.append((current_header, "".join(current_seq)))
                current_header = line
                current_seq = []
            elif line:
                current_seq.append(line)

        if current_header:
            sequences.append((current_header, "".join(current_seq)))

    return sequences


def merge_a3m_files(input_files: list[Path], output_path: Path):
    """Merge multiple A3M files, removing duplicates."""
    seen_seqs = set()
    all_sequences = []

    for input_file in input_files:
        sequences = parse_a3m(input_file)

        for header, seq in sequences:
            # Normalize for deduplication
            seq_key = seq.replace("-", "").replace(".", "").upper()

            if seq_key not in seen_seqs:
                seen_seqs.add(seq_key)
                all_sequences.append((header, seq))

    # Write merged output
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        for header, seq in all_sequences:
            f.write(f"{header}\n{seq}\n")

    logger.info(
        f"Merged {len(input_files)} files -> {len(all_sequences)} unique sequences"
    )
    logger.info(f"Output: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Merge multiple A3M files")
    parser.add_argument(
        "--inputs", nargs="+", type=Path, required=True, help="Input A3M files to merge"
    )
    parser.add_argument(
        "--output", type=Path, required=True, help="Output merged A3M file"
    )

    args = parser.parse_args()

    merge_a3m_files(args.inputs, args.output)


if __name__ == "__main__":
    main()
