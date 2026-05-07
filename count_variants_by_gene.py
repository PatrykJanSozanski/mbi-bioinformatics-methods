#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path


LOCAL_DEPS = Path(__file__).resolve().parent / "python_deps"
if LOCAL_DEPS.is_dir():
    sys.path.insert(0, str(LOCAL_DEPS))

import pandas as pd
import pyranges as pr


def read_refflat(path: Path) -> pd.DataFrame:
    columns = [
        "Gene",
        "Transcript",
        "Chromosome",
        "Strand",
        "Start",
        "End",
        "CdsStart",
        "CdsEnd",
        "ExonCount",
        "ExonStarts",
        "ExonEnds",
    ]
    df = pd.read_csv(path, sep="\t", header=None, names=columns)
    return df[["Chromosome", "Start", "End", "Strand", "Gene", "Transcript"]]


def read_vcf_variants(path: Path) -> pd.DataFrame:
    rows = []
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            rows.append(
                {
                    "Chromosome": chrom,
                    "Start": pos - 1,
                    "End": pos,
                    "VariantID": f"{chrom}:{pos}:{fields[3]}>{fields[4]}",
                }
            )
    return pd.DataFrame(rows, columns=["Chromosome", "Start", "End", "VariantID"])


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Count VCF variants overlapping genes from refFlat txStart-txEnd ranges."
    )
    parser.add_argument("vcf", type=Path)
    parser.add_argument("refflat", type=Path)
    parser.add_argument("-o", "--output", type=Path, default=None)
    parser.add_argument(
        "--include-zero",
        action="store_true",
        help="Include genes with no overlapping variants.",
    )
    args = parser.parse_args()

    genes_df = read_refflat(args.refflat)
    variants_df = read_vcf_variants(args.vcf)

    genes = pr.PyRanges(genes_df)
    variants = pr.PyRanges(variants_df)
    overlaps = variants.join(genes).df

    if overlaps.empty:
        counts = pd.DataFrame(columns=["Gene", "VariantCount"])
    else:
        unique_hits = overlaps[["Gene", "VariantID"]].drop_duplicates()
        counts = (
            unique_hits.groupby("Gene", as_index=False)
            .size()
            .rename(columns={"size": "VariantCount"})
        )

    if args.include_zero:
        all_genes = genes_df[["Gene"]].drop_duplicates()
        counts = all_genes.merge(counts, on="Gene", how="left").fillna(0)
        counts["VariantCount"] = counts["VariantCount"].astype(int)

    counts = counts.sort_values(["VariantCount", "Gene"], ascending=[False, True])

    if args.output:
        counts.to_csv(args.output, sep="\t", index=False)
    else:
        print(counts.to_csv(sep="\t", index=False), end="")


if __name__ == "__main__":
    main()

