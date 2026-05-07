#!/usr/bin/env python3

from argparse import ArgumentParser, ArgumentTypeError
from pathlib import Path

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


def fasta_path(value: str) -> Path:
    path = Path(value)
    if not path.is_file():
        raise ArgumentTypeError(f"{value} is not a valid path to a FASTA file")
    return path


def main() -> None:
    parser = ArgumentParser(
        description="Calculate GC-content for sequences from a FASTA file."
    )
    parser.add_argument("fasta", type=fasta_path, help="Path to FASTA file")
    args = parser.parse_args()

    total_gc_weight = 0.0
    total_length = 0
    sequence_count = 0

    for record in SeqIO.parse(args.fasta, "fasta"):
        sequence_count += 1
        seq_len = len(record.seq)
        gc_value = 100 * gc_fraction(record.seq)
        total_gc_weight += gc_fraction(record.seq) * seq_len
        total_length += seq_len

        print(f"Sequence ID: {record.id}")
        print(f"Length: {seq_len}")
        print(f"GC-content: {gc_value:.4f} %")

    if sequence_count == 0 or total_length == 0:
        raise SystemExit("No sequences found in the provided FASTA file.")

    weighted_gc = 100 * total_gc_weight / total_length
    print(f"Sequences processed: {sequence_count}")
    print(f"Total length: {total_length}")
    print(f"Weighted GC-content: {weighted_gc:.4f} %")


if __name__ == "__main__":
    main()
