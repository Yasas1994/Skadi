#!/usr/bin/env python3
"""
Compute all-vs-all pairwise scores and derive per-group thresholds.

This script processes an MMseqs2 all-vs-all search result (m8 file) to:
1. Compute per-pair ANI scores (tANI = ANI * symmetric coverage)
2. Label each pair by shared taxonomy ranks (species, genus, family, etc.)
3. For each taxonomic group and rank, fit score distributions and derive thresholds
4. Output thresholds as JSON for use in SKADI postprocessing

Usage:
    # Process existing m8 file
    python compute_pairwise_matrix.py /path/to/db /path/to/output \
        --m8 /path/to/all_vs_all.m8 --all

    # Or run MMseqs2 search first, then process
    python compute_pairwise_matrix.py /path/to/db /path/to/output \
        --run-search --all

Output:
    - pairs_labeled.parquet: All pairs with scores and taxonomy labels
    - thresholds.json: Per-group thresholds for each rank
    - threshold_report.txt: Human-readable summary
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import polars as pl
from scipy.stats import gaussian_kde
from tqdm import tqdm

from skadi.utils import ani_summary, axi_summary
from skadi.color_logger import logger


MMSEQS_HEADER = [
    "query", "target", "theader", "fident", "qlen", "tlen", "alnlen",
    "mismatch", "gapopen", "qstart", "qend", "tstart", "tend",
    "evalue", "bits", "taxid", "taxname", "taxlineage",
]

RANKS = ["species", "genus", "family", "order", "class", "phylum", "kingdom", "realm"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compute pairwise scores and derive per-group thresholds"
    )
    p.add_argument("db_dir", type=Path, help="Path to SKADI database directory")
    p.add_argument("outdir", type=Path, help="Output directory")
    p.add_argument("--m8", type=Path, default=None,
                   help="Path to existing MMseqs2 all-vs-all nucleotide m8 file")
    p.add_argument("--prot-m8", type=Path, default=None,
                   help="Path to protein search m8 file (for AAI threshold derivation)")
    p.add_argument("--prof-m8", type=Path, default=None,
                   help="Path to profile search m8 file (for API threshold derivation)")
    p.add_argument("--run-search", action="store_true",
                   help="Run MMseqs2 all-vs-all search (requires mmseqs in PATH)")
    p.add_argument("--compute-ani", action="store_true",
                   help="Compute ANI scores from m8")
    p.add_argument("--compute-aai", action="store_true",
                   help="Compute AAI scores from prot-m8")
    p.add_argument("--compute-api", action="store_true",
                   help="Compute API scores from prof-m8")
    p.add_argument("--label-pairs", action="store_true",
                   help="Label pairs by taxonomy")
    p.add_argument("--derive-thresholds", action="store_true",
                   help="Derive thresholds from labeled pairs")
    p.add_argument("--all", action="store_true",
                   help="Run all steps")
    p.add_argument("--threads", type=int, default=8,
                   help="Threads for MMseqs2 (default: 8)")
    p.add_argument("--min-pairs-per-group", type=int, default=10,
                   help="Minimum pairs required to derive threshold (default: 10)")
    p.add_argument("--threshold-method", type=str, default="youden",
                   choices=["youden", "percentile", "kde_crossing"],
                   help="Method for threshold computation (default: youden)")
    p.add_argument("--percentile", type=float, default=5.0,
                   help="Percentile for percentile method (default: 5)")
    p.add_argument("--score-type", type=str, default="tani",
                   choices=["tani", "ani", "aai", "api"],
                   help="Score type to use for threshold derivation (default: tani)")
    p.add_argument("--threshold-correction", type=float, default=1.0,
                   help="Scale all derived thresholds by this factor (default: 1.0)")
    p.add_argument("--max-threshold", type=float, default=None,
                   help="Cap thresholds at this maximum value (default: no cap)")
    p.add_argument("--use-reliable-only", action="store_true",
                   help="Only use thresholds marked reliable; fall back to global for others")
    return p.parse_args()


def run_mmseqs_all_vs_all(db_dir: Path, outdir: Path, threads: int = 8) -> Path:
    """Run MMseqs2 all-vs-all search on the genome database.

    Requires mmseqs to be in PATH and the database to have been indexed.
    """
    db_path = db_dir / "VMR_latest" / "genomes_fna"
    if not db_path.exists():
        # Try to find the database
        candidates = list((db_dir / "VMR_latest").glob("*"))
        logger.error("Database not found at %s. Candidates: %s", db_path, candidates)
        raise FileNotFoundError(f"MMseqs2 database not found: {db_path}")

    m8_out = outdir / "all_vs_all.m8"
    if m8_out.exists():
        logger.info("Using existing m8 file: %s", m8_out)
        return m8_out

    logger.info("Running MMseqs2 all-vs-all search...")
    cmd = [
        "mmseqs", "search",
        str(db_path), str(db_path),
        str(outdir / "result_db"), str(outdir / "tmp"),
        "--threads", str(threads),
        "-a", "--max-seqs", "1000",
        "--min-seq-id", "0.0", "-c", "0.0",
        "--cov-mode", "0", "--e", "10",
    ]
    subprocess.run(cmd, check=True)

    # Convert to m8
    cmd = [
        "mmseqs", "convertalis",
        str(db_path), str(db_path),
        str(outdir / "result_db"), str(m8_out),
        "--format-output", "query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage",
    ]
    subprocess.run(cmd, check=True)

    return m8_out


def compute_ani_scores(m8_path: Path, out_path: Path) -> pl.DataFrame:
    """Compute ANI scores from m8 file using ani_summary."""
    logger.info("Computing ANI scores from %s...", m8_path)
    df = ani_summary(
        infile=str(m8_path),
        all=True,
        header=MMSEQS_HEADER,
    )
    df.write_parquet(out_path)
    logger.info("Saved ANI scores to %s (%d pairs)", out_path, len(df))
    return df


def compute_aai_scores(
    prot_m8_path: Path,
    db_dir: Path,
    out_path: Path,
    chunksize: int = 1_000_000,
) -> pl.DataFrame:
    """Compute AAI scores from protein m8 file.

    Processes protein search results to compute tAAI per (query_genome, target_taxid).
    tAAI = AAI * qcov where qcov = unique_hit_proteins / total_query_proteins.

    Args:
        prot_m8_path: Path to protein search m8 file.
        db_dir: Database directory (for genome2protein mapping).
        out_path: Output parquet path.
        chunksize: Number of rows to process at a time.

    Returns:
        DataFrame with columns [seqid, taxid, aai, qcov, taai, hit_proteins, total_proteins].
    """
    logger.info("Computing AAI scores from %s...", prot_m8_path)

    # Load protein counts per genome
    g2p_path = db_dir / "VMR_latest" / "genome2protein"
    if not g2p_path.exists():
        logger.error("genome2protein file not found: %s", g2p_path)
        raise FileNotFoundError(g2p_path)

    g2p = pl.read_csv(g2p_path, separator="\t")
    protein_counts = g2p.group_by("genome_id").len().rename({"len": "total_proteins"})
    logger.info("Loaded protein counts for %d genomes", len(protein_counts))

    # Process m8 in chunks using pandas for chunked reading
    import pandas as pd

    header = [
        "query", "target", "theader", "fident", "qlen", "tlen", "alnlen",
        "mismatch", "gapopen", "qstart", "qend", "tstart", "tend",
        "evalue", "bits", "taxid", "taxname", "taxlineage",
    ]

    all_scores = []
    chunk_iter = pd.read_csv(
        prot_m8_path,
        sep="\t",
        header=None,
        names=header,
        chunksize=chunksize,
    )

    for i, pd_chunk in enumerate(chunk_iter):
        chunk = pl.from_pandas(pd_chunk)

        # Extract genome ID from query protein ID (format: genome_acc_ORFnum)
        chunk = chunk.with_columns(
            pl.col("query").str.extract(r"^(.+)_[0-9]+$", 1).alias("seqid")
        )

        # Join total proteins
        chunk = chunk.join(protein_counts, left_on="seqid", right_on="genome_id", how="left")

        # Compute AAI metrics per (seqid, taxid)
        scores = chunk.group_by(["seqid", "taxid"]).agg(
            ((pl.col("fident") * pl.col("alnlen")).sum() / pl.col("alnlen").sum())
            .round(3).alias("aai"),
            pl.col("query").n_unique().alias("hit_proteins"),
            pl.first("total_proteins").alias("total_proteins"),
        )

        # Compute qcov and taai
        scores = scores.with_columns(
            (pl.col("hit_proteins") / pl.col("total_proteins")).round(3).alias("qcov"),
        ).with_columns(
            (pl.col("aai") * pl.col("qcov")).round(3).alias("taai"),
        )

        all_scores.append(scores)
        if (i + 1) % 10 == 0:
            logger.info("  Processed %d chunks...", i + 1)

    if all_scores:
        df = pl.concat(all_scores)
    else:
        df = pl.DataFrame({
            "seqid": [], "taxid": [], "aai": [], "qcov": [], "taai": [],
            "hit_proteins": [], "total_proteins": [],
        })

    df.write_parquet(out_path)
    logger.info("Saved AAI scores to %s (%d pairs)", out_path, len(df))
    return df


def load_taxonomy(db_dir: Path) -> Tuple[Dict, Dict]:
    """Load taxonomy database and build name->taxid mapping.

    Returns:
        (taxdb, name2taxid) where taxdb is a taxopy.TaxDb object
    """
    from taxopy.core import TaxDb

    taxdb = TaxDb(
        nodes_dmp=f"{db_dir}/ictv-taxdump/nodes.dmp",
        names_dmp=f"{db_dir}/ictv-taxdump/names.dmp",
        merged_dmp=f"{db_dir}/ictv-taxdump/merged.dmp",
    )
    name2taxid = {v: k for k, v in taxdb.taxid2name.items()}
    return taxdb, name2taxid


def build_accession_taxonomy_map(
    db_dir: Path,
    taxdb,
    name2taxid: Dict,
) -> Dict[str, Dict[str, Optional[int]]]:
    """Build mapping from genome accession to rank taxids.

    Uses the ictv.cleaned.tsv file to map accessions to species names,
    then walks up the taxonomy tree to get taxids for each rank.
    """
    cleaned_path = db_dir / "ictv.cleaned.tsv"
    if not cleaned_path.exists():
        logger.warning("ictv.cleaned.tsv not found, trying ictv.tsv")
        cleaned_path = db_dir / "ictv.tsv"

    if not cleaned_path.exists():
        logger.error("No taxonomy file found in %s", db_dir)
        return {}

    df = pl.read_csv(cleaned_path, separator="\t")
    acc_col = "Virus GENBANK accession" if "Virus GENBANK accession" in df.columns else None
    species_col = "Species" if "Species" in df.columns else None

    if acc_col is None or species_col is None:
        logger.error("Required columns not found in %s", cleaned_path)
        return {}

    from taxopy import Taxon

    acc_to_ranks: Dict[str, Dict[str, Optional[int]]] = {}
    missing = 0

    for row in tqdm(df.iter_rows(named=True), total=len(df), desc="Building accession map"):
        acc = row[acc_col]
        species_name = row[species_col]
        if not acc or not species_name:
            continue

        sp_taxid = name2taxid.get(species_name)
        if not sp_taxid:
            missing += 1
            continue

        try:
            taxon = Taxon(sp_taxid, taxdb)
            rank_map = {taxdb.taxid2rank.get(t, "unknown"): t for t in taxon.taxid_lineage}
            acc_to_ranks[acc] = {r: rank_map.get(r) for r in RANKS}
        except Exception:
            pass

    logger.info("Mapped %d accessions (%d missing species names)", len(acc_to_ranks), missing)
    return acc_to_ranks


def label_pairs(
    pairs_df: pl.DataFrame,
    acc_to_ranks: Dict[str, Dict[str, Optional[int]]],
    out_path: Path,
) -> pl.DataFrame:
    """Label each pair by shared taxonomy ranks.

    For each pair (query, target), adds columns:
    - q_species_taxid, t_species_taxid, same_species, etc. for each rank
    """
    logger.info("Labeling pairs by taxonomy...")

    # Build lookup dataframe
    lookup_rows = []
    for acc, rank_map in acc_to_ranks.items():
        row = {"accession": acc}
        for r in RANKS:
            row[f"{r}_taxid"] = rank_map.get(r)
        lookup_rows.append(row)

    lookup = pl.DataFrame(lookup_rows)

    # Join query ranks
    pairs_q = pairs_df.join(
        lookup.rename(lambda c: f"q_{c}" if c != "accession" else "query"),
        on="query",
        how="left",
    )

    # Join target ranks
    pairs_labeled = pairs_q.join(
        lookup.rename(lambda c: f"t_{c}" if c != "accession" else "target"),
        on="target",
        how="left",
    )

    # Compute same-rank flags
    for rank in RANKS:
        pairs_labeled = pairs_labeled.with_columns(
            pl.when(
                pl.col(f"q_{rank}_taxid").is_not_null()
                & pl.col(f"t_{rank}_taxid").is_not_null()
                & (pl.col(f"q_{rank}_taxid") == pl.col(f"t_{rank}_taxid"))
            )
            .then(pl.lit(True))
            .otherwise(pl.lit(False))
            .alias(f"same_{rank}")
        )

    pairs_labeled.write_parquet(out_path)
    logger.info("Saved labeled pairs to %s", out_path)

    # Report statistics
    for rank in RANKS:
        same_count = pairs_labeled.filter(pl.col(f"same_{rank}")).height
        logger.info("  same_%s: %d (%.2f%%)", rank, same_count, 100 * same_count / len(pairs_labeled))

    return pairs_labeled


def fit_kde(scores: np.ndarray, bandwidth: Optional[float] = None) -> gaussian_kde:
    """Fit a Gaussian KDE to a score distribution."""
    if len(scores) < 2:
        return None
    # Filter out exact zeros for log-space stability if needed
    scores = scores[scores >= 0]
    if len(scores) < 2:
        return None
    try:
        kde = gaussian_kde(scores, bw_method=bandwidth or "scott")
        return kde
    except Exception as e:
        logger.warning("KDE fitting failed: %s", e)
        return None


def compute_threshold_youden(
    pos_scores: np.ndarray,
    neg_scores: np.ndarray,
    n_grid: int = 100,
) -> Tuple[float, float]:
    """Compute threshold maximizing Youden's J statistic.

    Returns (threshold, j_statistic).
    """
    if len(pos_scores) == 0 or len(neg_scores) == 0:
        return 0.0, 0.0

    all_scores = np.concatenate([pos_scores, neg_scores])
    min_s, max_s = all_scores.min(), all_scores.max()
    if min_s >= max_s:
        return min_s, 0.0

    thresholds = np.linspace(min_s, max_s, n_grid)
    best_j = -1.0
    best_t = thresholds[0]

    for t in thresholds:
        tp = np.sum(pos_scores >= t)
        fn = len(pos_scores) - tp
        tn = np.sum(neg_scores < t)
        fp = len(neg_scores) - tn

        sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
        j = sens + spec - 1.0

        if j > best_j:
            best_j = j
            best_t = t

    return float(best_t), float(best_j)


def compute_threshold_percentile(
    pos_scores: np.ndarray,
    percentile: float = 5.0,
) -> float:
    """Compute threshold as a percentile of positive scores."""
    if len(pos_scores) == 0:
        return 0.0
    return float(np.percentile(pos_scores, percentile))


def compute_threshold_kde_crossing(
    pos_scores: np.ndarray,
    neg_scores: np.ndarray,
    n_grid: int = 500,
) -> Tuple[float, float]:
    """Compute threshold as the crossing point of KDEs.

    Finds the highest score where pos_kde > neg_kde.
    """
    if len(pos_scores) < 3 or len(neg_scores) < 3:
        return 0.0, 0.0

    kde_pos = fit_kde(pos_scores)
    kde_neg = fit_kde(neg_scores)

    if kde_pos is None or kde_neg is None:
        return 0.0, 0.0

    all_scores = np.concatenate([pos_scores, neg_scores])
    min_s, max_s = all_scores.min(), all_scores.max()
    if min_s >= max_s:
        return min_s, 0.0

    x = np.linspace(max(min_s, 0.001), max_s, n_grid)
    pos_pdf = kde_pos(x)
    neg_pdf = kde_neg(x)

    # Find crossing from right to left (highest x where pos > neg)
    crossing = None
    for i in range(len(x) - 1, -1, -1):
        if pos_pdf[i] > neg_pdf[i]:
            crossing = x[i]
            break

    if crossing is None:
        crossing = float(np.percentile(pos_scores, 5))

    return crossing, 0.0


def derive_group_thresholds(
    pairs_df: pl.DataFrame,
    group_rank: str,
    target_rank: str,
    score_col: str,
    method: str = "youden",
    percentile: float = 5.0,
    min_pairs: int = 10,
    taxdb=None,
) -> Dict:
    """Derive thresholds for a specific group rank and target rank.

    For example, group_rank="family", target_rank="genus" means:
    - For each family, find the threshold that best separates same-genus from diff-genus pairs.

    Args:
        pairs_df: Labeled pairs dataframe
        group_rank: Rank to group by (e.g., "family")
        target_rank: Rank to derive threshold for (e.g., "genus")
        score_col: Column name for scores (e.g., "tani")
        method: Threshold computation method
        min_pairs: Minimum pairs required per group

    Returns:
        Dict mapping group taxid -> threshold info
    """
    results = {}

    # Get unique groups
    group_col = f"q_{group_rank}_taxid"
    groups = pairs_df.filter(pl.col(group_col).is_not_null())[group_col].unique().to_list()

    for group_taxid in tqdm(groups, desc=f"{group_rank}->{target_rank}"):
        group_pairs = pairs_df.filter(pl.col(group_col) == group_taxid)

        # Positives: pairs that share the target rank
        # For genus within family: same_genus = True
        pos = group_pairs.filter(pl.col(f"same_{target_rank}"))

        # Negatives: pairs in the same group but different target rank
        # For genus within family: same_family but not same_genus
        if target_rank == group_rank:
            # For family threshold within family: compare same_family vs diff_family
            # But we need negatives from outside the group
            neg = pairs_df.filter(
                (pl.col(group_col) != group_taxid) | pl.col(group_col).is_null()
            )
        else:
            neg = group_pairs.filter(~pl.col(f"same_{target_rank}"))

        if pos.height < min_pairs:
            continue

        pos_scores = pos[score_col].to_numpy()
        neg_scores = neg[score_col].to_numpy() if neg.height > 0 else np.array([])

        if method == "youden":
            thresh, j = compute_threshold_youden(pos_scores, neg_scores)
        elif method == "percentile":
            thresh = compute_threshold_percentile(pos_scores, percentile)
            j = 0.0
        elif method == "kde_crossing":
            thresh, j = compute_threshold_kde_crossing(pos_scores, neg_scores)
        else:
            thresh, j = 0.0, 0.0

        # Compute statistics
        sens = np.sum(pos_scores >= thresh) / len(pos_scores) if len(pos_scores) > 0 else 0.0
        spec = np.sum(neg_scores < thresh) / len(neg_scores) if len(neg_scores) > 0 else 1.0

        group_name = taxdb.taxid2name.get(group_taxid, "Unknown") if taxdb else "Unknown"

        # Apply sanity checks and fallback logic
        final_thresh = float(thresh)
        # If threshold is unrealistically high (>0.95) or low (<0.001),
        # mark it as potentially unreliable
        reliable = True
        if final_thresh > 0.95 and len(pos_scores) < 50:
            reliable = False
        if final_thresh < 0.001:
            reliable = False

        results[int(group_taxid)] = {
            "name": group_name,
            "threshold": round(final_thresh, 4),
            "sensitivity": round(float(sens), 4),
            "specificity": round(float(spec), 4),
            "j_statistic": round(float(j), 4),
            "n_positive": int(len(pos_scores)),
            "n_negative": int(len(neg_scores)),
            "pos_median": round(float(np.median(pos_scores)), 4),
            "neg_median": round(float(np.median(neg_scores)), 4) if len(neg_scores) > 0 else None,
            "reliable": reliable,
        }

    return results


def derive_thresholds_for_method(
    pairs_df: pl.DataFrame,
    score_col: str,
    target_ranks: List[str],
    method: str = "youden",
    percentile: float = 5.0,
    min_pairs: int = 10,
    taxdb=None,
) -> Dict:
    """Derive thresholds for specific target ranks from a score column.

    Args:
        pairs_df: Labeled pairs dataframe.
        score_col: Score column to use (e.g., "tani", "taai", "tapi").
        target_ranks: List of ranks to derive thresholds for.
        method: Threshold computation method.
        percentile: Percentile for percentile method.
        min_pairs: Minimum pairs per group.
        taxdb: TaxDb object for name lookups.

    Returns:
        Dict with structure: target_rank -> group_rank -> group_taxid -> threshold_info
    """
    results = {}

    # Map target rank to appropriate group rank
    rank_to_group = {
        "species": "family",
        "genus": "family",
        "family": "order",
        "order": "class",
        "class": "phylum",
        "phylum": "kingdom",
        "kingdom": "realm",
    }

    for target_rank in target_ranks:
        group_rank = rank_to_group.get(target_rank)
        if not group_rank:
            continue

        logger.info("Deriving %s thresholds per %s (%s)...", target_rank, group_rank, score_col)
        group_results = derive_group_thresholds(
            pairs_df,
            group_rank,
            target_rank,
            score_col,
            method=method,
            percentile=percentile,
            min_pairs=min_pairs,
            taxdb=taxdb,
        )

        if target_rank not in results:
            results[target_rank] = {}
        results[target_rank][group_rank] = group_results
        logger.info("  Derived %d thresholds", len(group_results))

    # Compute global thresholds for target ranks
    results["global"] = {}
    for target_rank in target_ranks:
        pos = pairs_df.filter(pl.col(f"same_{target_rank}"))
        neg = pairs_df.filter(~pl.col(f"same_{target_rank}"))

        if pos.height < min_pairs:
            continue

        pos_scores = pos[score_col].to_numpy()
        neg_scores = neg[score_col].to_numpy() if neg.height > 0 else np.array([])

        if method == "youden":
            thresh, j = compute_threshold_youden(pos_scores, neg_scores)
        elif method == "percentile":
            thresh = compute_threshold_percentile(pos_scores, percentile)
            j = 0.0
        elif method == "kde_crossing":
            thresh, j = compute_threshold_kde_crossing(pos_scores, neg_scores)
        else:
            thresh, j = 0.0, 0.0

        sens = np.sum(pos_scores >= thresh) / len(pos_scores) if len(pos_scores) > 0 else 0.0
        spec = np.sum(neg_scores < thresh) / len(neg_scores) if len(neg_scores) > 0 else 1.0

        results["global"][target_rank] = {
            "threshold": round(float(thresh), 4),
            "sensitivity": round(float(sens), 4),
            "specificity": round(float(spec), 4),
            "j_statistic": round(float(j), 4),
            "n_positive": int(len(pos_scores)),
            "n_negative": int(len(neg_scores)),
        }

    return results


def apply_threshold_adjustments(
    thresholds: Dict,
    correction: float = 1.0,
    max_threshold: Optional[float] = None,
    use_reliable_only: bool = False,
) -> Dict:
    """Apply correction factor, max cap, and reliability filtering to thresholds.

    Args:
        thresholds: Nested dict from derive_all_thresholds.
        correction: Multiply all thresholds by this factor.
        max_threshold: Cap thresholds at this value.
        use_reliable_only: If True, set unreliable per-group thresholds to None
                           (so downstream uses global fallback).

    Returns:
        Modified thresholds dict.
    """
    adjusted = {}

    for key, value in thresholds.items():
        if key == "global":
            adjusted[key] = {}
            for rank, info in value.items():
                new_info = dict(info)
                t = new_info["threshold"] * correction
                if max_threshold is not None:
                    t = min(t, max_threshold)
                new_info["threshold"] = round(float(t), 4)
                new_info["original_threshold"] = round(float(info["threshold"]), 4)
                adjusted[key][rank] = new_info
        elif isinstance(value, dict):
            adjusted[key] = {}
            for group_rank, groups in value.items():
                adjusted[key][group_rank] = {}
                for taxid, info in groups.items():
                    new_info = dict(info)
                    t = new_info["threshold"] * correction
                    if max_threshold is not None:
                        t = min(t, max_threshold)

                    # If use_reliable_only and not reliable, set to global fallback
                    if use_reliable_only and not new_info.get("reliable", True):
                        # Will be handled by downstream fallback logic
                        new_info["threshold"] = None
                        new_info["disabled"] = True
                    else:
                        new_info["threshold"] = round(float(t), 4)

                    new_info["original_threshold"] = round(float(info["threshold"]), 4)
                    adjusted[key][group_rank][taxid] = new_info
        else:
            adjusted[key] = value

    return adjusted


def write_threshold_report(thresholds: Dict, out_path: Path) -> None:
    """Write a human-readable threshold report."""
    with open(out_path, "w") as f:
        f.write("# SKADI Per-Group Threshold Report\n\n")

        # Global thresholds
        f.write("## Global Thresholds\n\n")
        f.write(f"{'Rank':<12} {'Threshold':>10} {'Sensitivity':>12} {'Specificity':>12} {'N+':>8} {'N-':>8}\n")
        f.write("-" * 60 + "\n")
        for rank, info in thresholds.get("global", {}).items():
            f.write(f"{rank:<12} {info['threshold']:>10.4f} {info['sensitivity']:>12.4f} "
                    f"{info['specificity']:>12.4f} {info['n_positive']:>8} {info['n_negative']:>8}\n")

        # Per-group thresholds
        for target_rank in ["species", "genus", "family", "order", "class", "phylum", "kingdom"]:
            if target_rank not in thresholds:
                continue

            for group_rank, groups in thresholds[target_rank].items():
                if not groups:
                    continue

                f.write(f"\n## {target_rank.capitalize()} Thresholds per {group_rank.capitalize()}\n\n")
                f.write(f"{'Group':<30} {'Threshold':>10} {'Sens':>8} {'Spec':>8} {'N+':>8} {'N-':>8}\n")
                f.write("-" * 80 + "\n")

                # Sort by threshold
                sorted_groups = sorted(groups.items(), key=lambda x: x[1]["threshold"])
                for taxid, info in sorted_groups:
                    name = info.get("name", "Unknown")[:28]
                    f.write(f"{name:<30} {info['threshold']:>10.4f} {info['sensitivity']:>8.4f} "
                            f"{info['specificity']:>8.4f} {info['n_positive']:>8} {info['n_negative']:>8}\n")

    logger.info("Threshold report saved to %s", out_path)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.all:
        args.compute_ani = True
        args.label_pairs = True
        args.derive_thresholds = True

    # Step 1: Get or run MMseqs2 search
    m8_path = args.m8
    if args.run_search:
        m8_path = run_mmseqs_all_vs_all(args.db_dir, args.outdir, args.threads)

    # Check if we can skip m8 requirement because labeled pairs already exist
    labeled_path = args.outdir / "pairs_labeled.parquet"
    ani_path = args.outdir / "pairs_ani.parquet"

    if m8_path is None:
        # Check if there's an existing m8 in the output dir
        existing = list(args.outdir.glob("*.m8"))
        if existing:
            m8_path = existing[0]
            logger.info("Using existing m8 file: %s", m8_path)
        elif labeled_path.exists() and not (args.compute_ani or args.run_search):
            # Can skip m8 if we already have labeled pairs and don't need to recompute
            logger.info("Using existing labeled pairs: %s", labeled_path)
        elif ani_path.exists() and not (args.compute_ani or args.run_search):
            # Can skip m8 if we already have ANI scores and just need to label
            logger.info("Using existing ANI scores: %s", ani_path)
        else:
            logger.error("No m8 file provided. Use --m8 or --run-search.")
            sys.exit(1)

    # Step 2: Compute ANI scores
    ani_path = args.outdir / "pairs_ani.parquet"
    if args.compute_ani or args.all:
        if not ani_path.exists():
            compute_ani_scores(m8_path, ani_path)
        else:
            logger.info("Using existing ANI scores: %s", ani_path)

    # Step 3: Label pairs by taxonomy
    labeled_path = args.outdir / "pairs_labeled.parquet"
    if args.label_pairs or args.all:
        if not labeled_path.exists():
            taxdb, name2taxid = load_taxonomy(args.db_dir)
            acc_to_ranks = build_accession_taxonomy_map(args.db_dir, taxdb, name2taxid)

            pairs_df = pl.read_parquet(ani_path)
            label_pairs(pairs_df, acc_to_ranks, labeled_path)
        else:
            logger.info("Using existing labeled pairs: %s", labeled_path)

    # Step 4: Compute AAI scores if protein m8 provided
    aai_path = args.outdir / "pairs_aai.parquet"
    if args.compute_aai and args.prot_m8:
        if not aai_path.exists():
            compute_aai_scores(args.prot_m8, args.db_dir, aai_path)
        else:
            logger.info("Using existing AAI scores: %s", aai_path)

    # Step 5: Label AAI pairs by taxonomy
    aai_labeled_path = args.outdir / "pairs_aai_labeled.parquet"
    if args.compute_aai and args.prot_m8 and not aai_labeled_path.exists():
        taxdb, name2taxid = load_taxonomy(args.db_dir)
        acc_to_ranks = build_accession_taxonomy_map(args.db_dir, taxdb, name2taxid)
        aai_df = pl.read_parquet(aai_path)
        # Rename seqid to query for consistency with label_pairs
        aai_df = aai_df.rename({"seqid": "query"})
        label_pairs(aai_df, acc_to_ranks, aai_labeled_path)

    # Step 6: Derive thresholds
    if args.derive_thresholds or args.all:
        taxdb, _ = load_taxonomy(args.db_dir)
        all_thresholds = {}

        # ANI thresholds: species + genus
        if ani_path.exists():
            logger.info("Deriving ANI thresholds (species, genus)...")
            ani_df = pl.read_parquet(labeled_path)
            ani_thresholds = derive_thresholds_for_method(
                ani_df,
                score_col="tani",
                target_ranks=["species", "genus"],
                method=args.threshold_method,
                percentile=args.percentile,
                min_pairs=args.min_pairs_per_group,
                taxdb=taxdb,
            )
            all_thresholds["ani"] = ani_thresholds

        # AAI thresholds: family + order + class + phylum + kingdom
        if aai_labeled_path.exists():
            logger.info("Deriving AAI thresholds (family, order, class, phylum, kingdom)...")
            aai_df = pl.read_parquet(aai_labeled_path)
            aai_thresholds = derive_thresholds_for_method(
                aai_df,
                score_col="taai",
                target_ranks=["family", "order", "class", "phylum", "kingdom"],
                method=args.threshold_method,
                percentile=args.percentile,
                min_pairs=args.min_pairs_per_group,
                taxdb=taxdb,
            )
            all_thresholds["aai"] = aai_thresholds

        # API thresholds: family + order + class + phylum + kingdom
        # TODO: Implement when profile search data is available
        if args.prof_m8:
            logger.warning("API threshold derivation not yet implemented")

        # Apply threshold correction and caps
        if args.threshold_correction != 1.0 or args.max_threshold is not None:
            logger.info("Applying threshold correction (factor=%.2f, max=%s)...",
                        args.threshold_correction, args.max_threshold)
            for method_key in all_thresholds:
                all_thresholds[method_key] = apply_threshold_adjustments(
                    all_thresholds[method_key],
                    correction=args.threshold_correction,
                    max_threshold=args.max_threshold,
                    use_reliable_only=args.use_reliable_only,
                )

        # Save thresholds
        thresholds_path = args.outdir / "thresholds.json"
        with open(thresholds_path, "w") as f:
            json.dump(all_thresholds, f, indent=2)
        logger.info("Thresholds saved to %s", thresholds_path)

        # Write report
        write_threshold_report(all_thresholds, args.outdir / "threshold_report.txt")


if __name__ == "__main__":
    main()
