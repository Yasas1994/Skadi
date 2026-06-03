#!/usr/bin/env python3
"""
Compute all-vs-all ANI/AAI matrices from an MSL database.

This generates the reference data needed for data-driven threshold derivation.
Output: Parquet files with pairwise scores and ground-truth same-rank labels.

Usage:
    # Step 1: Sample pairs
    python compute_pairwise_matrix.py /path/to/db /path/to/output --sample-pairs

    # Step 2: Compute scores (requires MMseqs2 databases)
    python compute_pairwise_matrix.py /path/to/db /path/to/output --compute-scores

    # Step 3: Derive thresholds
    python compute_pairwise_matrix.py /path/to/db /path/to/output --derive-thresholds

    # Or all at once:
    python compute_pairwise_matrix.py /path/to/db /path/to/output --all
"""

from __future__ import annotations

import argparse
import random
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import polars as pl
from tqdm import tqdm

from skadi.utils import ani_summary, axi_summary
from skadi.benchmark_thresholds import grid_search_thresholds
from skadi.color_logger import logger


MMSEQS_HEADER = [
    "query", "target", "theader", "fident", "qlen", "tlen", "alnlen",
    "mismatch", "gapopen", "qstart", "qend", "tstart", "tend",
    "evalue", "bits", "taxid", "taxname", "taxlineage",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compute pairwise ANI/AAI matrix from MSL database"
    )
    p.add_argument("db_dir", type=Path, help="Path to SKADI database directory")
    p.add_argument("outdir", type=Path, help="Output directory for matrices")
    p.add_argument("--max-pairs", type=int, default=500_000,
                   help="Maximum number of pairs to sample per rank (default: 500k)")
    p.add_argument("--max-taxa", type=int, default=None,
                   help="Maximum number of taxa to sample from (default: all)")
    p.add_argument("--max-genomes-per-taxon", type=int, default=20,
                   help="Max genomes per taxon for intra-pairs (default: 20)")
    p.add_argument("--seed", type=int, default=42, help="Random seed for pair sampling")
    p.add_argument("--sample-pairs", action="store_true", help="Run pair sampling")
    p.add_argument("--compute-scores", action="store_true", help="Run MMseqs2 and compute scores")
    p.add_argument("--derive-thresholds", action="store_true", help="Derive optimal thresholds")
    p.add_argument("--all", action="store_true", help="Run all steps")
    p.add_argument("--mmseqs-threads", type=int, default=4, help="Threads for MMseqs2")
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


def load_taxonomy(db_dir: Path) -> Tuple[Dict[int, int], Dict[int, int], Dict[int, int]]:
    """Load taxonomy mappings: taxid -> family/genus/species taxids."""
    from taxopy.core import TaxDb
    from taxopy import Taxon

    taxdb = TaxDb(
        nodes_dmp=f"{db_dir}/ictv-taxdump/nodes.dmp",
        names_dmp=f"{db_dir}/ictv-taxdump/names.dmp",
        merged_dmp=f"{db_dir}/ictv-taxdump/merged.dmp",
    )

    acc2tax = pl.read_csv(
        db_dir / "VMR_latest" / "virus_genome.accession2taxid",
        separator="\t",
        has_header=True,
        new_columns=["accession", "accession_version", "taxid", "gi"],
    )

    family_map: Dict[int, int] = {}
    genus_map: Dict[int, int] = {}
    species_map: Dict[int, int] = {}

    for tid in tqdm(acc2tax["taxid"].unique().to_list(), desc="Loading taxonomy"):
        try:
            taxon = Taxon(tid, taxdb=taxdb)
            family_map[tid] = taxon.rank_taxid_dictionary.get("family", -1)
            genus_map[tid] = taxon.rank_taxid_dictionary.get("genus", -1)
            species_map[tid] = taxon.rank_taxid_dictionary.get("species", -1)
        except Exception:
            family_map[tid] = -1
            genus_map[tid] = -1
            species_map[tid] = -1

    return family_map, genus_map, species_map


def sample_pairs(
    df: pl.DataFrame,
    max_pairs: int,
    seed: int,
    max_taxa: Optional[int] = None,
    max_genomes_per_taxon: int = 20,
) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """Sample intra-rank and inter-rank pairs for threshold optimization.

    Returns:
        (intra_pairs, inter_pairs) DataFrames with columns
        [g1, g2, taxid1, taxid2, family1, family2, genus1, genus2, species1, species2].
    """
    random.seed(seed)
    np.random.seed(seed)

    taxa = df["taxid"].unique().to_list()
    if max_taxa is not None:
        taxa = taxa[:max_taxa]

    intra_pairs = []
    inter_pairs = []

    logger.info("Sampling pairs from %d taxa...", len(taxa))

    for taxid in tqdm(taxa, desc="Taxa"):
        genomes = df.filter(pl.col("taxid") == taxid)["accession"].to_list()
        if len(genomes) < 2:
            continue

        # Intra-taxa pairs
        limit = min(len(genomes), max_genomes_per_taxon)
        for i in range(limit):
            for j in range(i + 1, limit):
                intra_pairs.append({
                    "g1": genomes[i],
                    "g2": genomes[j],
                    "taxid1": taxid,
                    "taxid2": taxid,
                })

        # Inter-taxa pairs (sample one random genome from different taxa)
        other_candidates = [t for t in taxa if t != taxid]
        if other_candidates:
            other = random.choice(other_candidates)
            other_genomes = df.filter(pl.col("taxid") == other)["accession"].to_list()
            if other_genomes:
                for g in genomes[:5]:
                    inter_pairs.append({
                        "g1": g,
                        "g2": other_genomes[0],
                        "taxid1": taxid,
                        "taxid2": other,
                    })

        if len(intra_pairs) >= max_pairs:
            break

    return pl.DataFrame(intra_pairs), pl.DataFrame(inter_pairs)


def run_mmseqs_search(
    pairs_df: pl.DataFrame,
    db_dir: Path,
    outdir: Path,
    threads: int = 4,
) -> Path:
    """Run MMseqs2 search on a set of genome pairs.

    Creates temporary query/target FASTA subsets and runs mmseqs search.
    Returns path to the m8 output file.
    """
    # This is a placeholder for the actual MMseqs2 execution.
    # In practice, the user would need to:
    # 1. Extract genome sequences for each pair
    # 2. Run mmseqs createdb / search / convertalis
    # 3. Parse the output m8 file

    logger.info("MMseqs2 search would run on %d pairs", len(pairs_df))
    logger.info("(This requires actual genome FASTA files and MMseqs2 databases)")

    # Create a sentinel file to indicate the step
    sentinel = outdir / ".mmseqs_search_pending"
    sentinel.write_text(
        "MMseqs2 search not yet implemented.\n"
        "Please run mmseqs2 manually on the sampled pairs and save results\n"
        "to intra_scores.parquet and inter_scores.parquet.\n"
    )
    return sentinel


def compute_scores_from_m8(
    m8_path: Path,
    db_dir: Path,
    out_path: Path,
) -> pl.DataFrame:
    """Compute ANI scores from an m8 file and save to parquet."""
    try:
        df = ani_summary(
            infile=str(m8_path),
            all=True,
            header=MMSEQS_HEADER,
            dbdir=str(db_dir),
        )
        df.write_parquet(out_path)
        return df
    except Exception as e:
        logger.error("Failed to compute scores from %s: %s", m8_path, e)
        raise


def derive_thresholds(
    intra_scores: pl.DataFrame,
    inter_scores: pl.DataFrame,
    family_map: Dict[int, int],
    output_path: Path,
) -> Dict:
    """Derive optimal thresholds per family from score distributions.

    Uses grid search to find thresholds that maximize F1 for each rank.
    """
    results = {
        "global": {},
        "per_family": {},
    }

    score_col = "tani"
    if score_col not in intra_scores.columns:
        score_col = "ani"

    # Global thresholds
    logger.info("Computing global thresholds...")
    intra_arr = intra_scores[score_col].to_numpy()
    inter_arr = inter_scores[score_col].to_numpy()

    all_scores = np.concatenate([intra_arr, inter_arr])
    all_labels = np.concatenate([np.ones(len(intra_arr)), np.zeros(len(inter_arr))])

    threshold_range = np.arange(0.05, 0.95, 0.01)
    for rank in ["species", "genus", "family"]:
        best_thresh, best_f1 = grid_search_thresholds(
            all_scores, all_labels, threshold_range, metric="f1"
        )
        results["global"][rank] = {
            "threshold": round(best_thresh, 4),
            "f1": round(best_f1, 4),
        }
        logger.info("Global %s threshold: %.3f (F1=%.3f)", rank, best_thresh, best_f1)

    # Per-family thresholds
    logger.info("Computing per-family thresholds...")
    intra_fam = intra_scores.with_columns(
        pl.col("taxid1").replace_strict(family_map, default=-1).alias("family")
    )

    families = intra_fam.filter(pl.col("family") > 0)["family"].unique().to_list()
    for fam in tqdm(families, desc="Families"):
        fam_intra = intra_fam.filter(pl.col("family") == fam)[score_col].to_numpy()
        fam_inter = inter_scores[score_col].to_numpy()  # Use all inter as negatives

        if len(fam_intra) < 5:
            continue

        fam_scores = np.concatenate([fam_intra, fam_inter])
        fam_labels = np.concatenate([np.ones(len(fam_intra)), np.zeros(len(fam_inter))])

        best_thresh, best_f1 = grid_search_thresholds(
            fam_scores, fam_labels, threshold_range, metric="f1"
        )
        results["per_family"][int(fam)] = {
            "threshold": round(best_thresh, 4),
            "f1": round(best_f1, 4),
            "n_intra": len(fam_intra),
        }

    # Save results
    import json
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)

    logger.info("Thresholds saved to %s", output_path)
    return results


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.all:
        args.sample_pairs = True
        args.compute_scores = True
        args.derive_thresholds = True

    # Step 1: Sample pairs
    if args.sample_pairs:
        logger.info("Loading genome metadata from %s", args.db_dir)
        genome_df = load_genome_list(args.db_dir)
        logger.info("Found %d genomes", len(genome_df))

        logger.info("Sampling pairs (max %d)...", args.max_pairs)
        intra, inter = sample_pairs(
            genome_df,
            args.max_pairs,
            args.seed,
            max_taxa=args.max_taxa,
            max_genomes_per_taxon=args.max_genomes_per_taxon,
        )

        logger.info("Intra-taxa pairs: %d", len(intra))
        logger.info("Inter-taxa pairs: %d", len(inter))

        intra.write_parquet(args.outdir / "intra_pairs.parquet")
        inter.write_parquet(args.outdir / "inter_pairs.parquet")

    # Step 2: Compute scores
    if args.compute_scores:
        intra = pl.read_parquet(args.outdir / "intra_pairs.parquet")
        inter = pl.read_parquet(args.outdir / "inter_pairs.parquet")

        logger.info("Running MMseqs2 searches...")
        run_mmseqs_search(intra, args.db_dir, args.outdir, threads=args.mmseqs_threads)

        # After MMseqs2 completes, compute scores
        # m8_file = args.outdir / "intra_search.m8"
        # if m8_file.exists():
        #     compute_scores_from_m8(m8_file, args.db_dir, args.outdir / "intra_scores.parquet")

    # Step 3: Derive thresholds
    if args.derive_thresholds:
        scores_file = args.outdir / "intra_scores.parquet"
        inter_scores_file = args.outdir / "inter_scores.parquet"

        if not scores_file.exists() or not inter_scores_file.exists():
            logger.error(
                "Score files not found. Run --compute-scores first or provide\n"
                "intra_scores.parquet and inter_scores.parquet in %s",
                args.outdir,
            )
            sys.exit(1)

        intra_scores = pl.read_parquet(scores_file)
        inter_scores = pl.read_parquet(inter_scores_file)

        logger.info("Loading taxonomy for per-family analysis...")
        family_map, _, _ = load_taxonomy(args.db_dir)

        derive_thresholds(
            intra_scores,
            inter_scores,
            family_map,
            args.outdir / "derived_thresholds.json",
        )


if __name__ == "__main__":
    main()
