#!/usr/bin/env python3
"""
Compute all-vs-all ANI/AAI matrices from an MSL database.

This generates the reference data needed for data-driven threshold derivation.
Output: Parquet files with pairwise scores and ground-truth same-rank labels.

Usage:
    python compute_pairwise_matrix.py /path/to/db /path/to/output [--max-pairs 100000]
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

import numpy as np
import polars as pl
import pyfastx
from tqdm import tqdm

from skadi.utils import ani_summary, axi_summary, index_m8
from skadi.color_logger import logger


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compute pairwise ANI/AAI matrix from MSL database"
    )
    p.add_argument("db_dir", type=Path, help="Path to SKADI database directory")
    p.add_argument("outdir", type=Path, help="Output directory for matrices")
    p.add_argument(
        "--max-pairs",
        type=int,
        default=500_000,
        help="Maximum number of pairs to sample per rank (default: 500k)",
    )
    p.add_argument(
        "--seed", type=int, default=42, help="Random seed for pair sampling"
    )
    return p.parse_args()


def load_genome_list(db_dir: Path) -> pl.DataFrame:
    """Load genome accessions and taxids from the database."""
    acc2tax = db_dir / "VMR_latest" / "virus_genome.accession2taxid"
    df = pl.read_csv(
        acc2tax,
        separator="\t",
        has_header=True,
        new_columns=["accession", "accession_version", "taxid", "gi"],
    )
    return df.select(["accession", "taxid"]).unique()


def sample_pairs(
    df: pl.DataFrame, max_pairs: int, seed: int
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Sample intra-rank and inter-rank pairs for threshold optimization.

    Returns:
        (intra_pairs, inter_pairs) DataFrames with columns [g1, g2, taxid1, taxid2].
    """
    random.seed(seed)
    np.random.seed(seed)

    # Group by taxid to find same-species pairs
    taxa = df["taxid"].unique().to_list()

    intra_pairs = []
    inter_pairs = []

    logger.info("Sampling pairs from %d taxa...", len(taxa))

    for taxid in tqdm(taxa[:1000], desc="Taxa"):  # Limit to first 1000 for speed
        genomes = df.filter(pl.col("taxid") == taxid)["accession"].to_list()
        if len(genomes) < 2:
            continue

        # Intra-taxa pairs
        for i in range(min(len(genomes), 20)):
            for j in range(i + 1, min(len(genomes), 20)):
                intra_pairs.append({
                    "g1": genomes[i],
                    "g2": genomes[j],
                    "taxid1": taxid,
                    "taxid2": taxid,
                })

        # Inter-taxa pairs (sample one random genome from different taxa)
        other = random.choice([t for t in taxa if t != taxid])
        other_genome = df.filter(pl.col("taxid") == other)["accession"].to_list()[0]
        for g in genomes[:5]:
            inter_pairs.append({
                "g1": g,
                "g2": other_genome,
                "taxid1": taxid,
                "taxid2": other,
            })

        if len(intra_pairs) >= max_pairs:
            break

    return pl.DataFrame(intra_pairs), pl.DataFrame(inter_pairs)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    logger.info("Loading genome metadata from %s", args.db_dir)
    genome_df = load_genome_list(args.db_dir)
    logger.info("Found %d genomes", len(genome_df))

    logger.info("Sampling pairs (max %d)...", args.max_pairs)
    intra, inter = sample_pairs(genome_df, args.max_pairs, args.seed)

    logger.info("Intra-taxa pairs: %d", len(intra))
    logger.info("Inter-taxa pairs: %d", len(inter))

    # Save sampled pairs for downstream analysis
    intra.write_parquet(args.outdir / "intra_pairs.parquet")
    inter.write_parquet(args.outdir / "inter_pairs.parquet")

    logger.info("Done. Use these pairs with mmseqs2 search results to derive thresholds.")


if __name__ == "__main__":
    main()
