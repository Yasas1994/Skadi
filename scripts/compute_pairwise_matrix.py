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
                   help="Path to existing MMseqs2 all-vs-all m8 file")
    p.add_argument("--run-search", action="store_true",
                   help="Run MMseqs2 all-vs-all search (requires mmseqs in PATH)")
    p.add_argument("--compute-ani", action="store_true",
                   help="Compute ANI scores from m8")
    p.add_argument("--compute-aai", action="store_true",
                   help="Compute AAI scores (requires protein search)")
    p.add_argument("--compute-api", action="store_true",
                   help="Compute API scores (requires protein search)")
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


def derive_all_thresholds(
    pairs_df: pl.DataFrame,
    score_col: str,
    method: str = "youden",
    percentile: float = 5.0,
    min_pairs: int = 10,
    taxdb=None,
) -> Dict:
    """Derive thresholds for all group/rank combinations.

    Returns nested dict: group_rank -> target_rank -> group_taxid -> threshold_info
    """
    results = {}

    # Define which group ranks to use for each target rank
    # For species: group by family (or use global)
    # For genus: group by family
    # For family: group by order (or class if order missing)
    # etc.

    group_target_pairs = [
        ("family", "species"),
        ("family", "genus"),
        ("order", "family"),
        ("class", "order"),
        ("phylum", "class"),
        ("kingdom", "phylum"),
        ("realm", "kingdom"),
    ]

    for group_rank, target_rank in group_target_pairs:
        logger.info("Deriving %s thresholds per %s...", target_rank, group_rank)
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

    # Also compute global thresholds
    logger.info("Computing global thresholds...")
    results["global"] = {}
    for target_rank in ["species", "genus", "family", "order", "class", "phylum", "kingdom"]:
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

    # Step 4: Derive thresholds
    if args.derive_thresholds or args.all:
        pairs_df = pl.read_parquet(labeled_path)
        taxdb, _ = load_taxonomy(args.db_dir)

        logger.info("Deriving thresholds using %s method...", args.threshold_method)
        thresholds = derive_all_thresholds(
            pairs_df,
            score_col=args.score_type,
            method=args.threshold_method,
            percentile=args.percentile,
            min_pairs=args.min_pairs_per_group,
            taxdb=taxdb,
        )

        # Save thresholds
        thresholds_path = args.outdir / "thresholds.json"
        with open(thresholds_path, "w") as f:
            json.dump(thresholds, f, indent=2)
        logger.info("Thresholds saved to %s", thresholds_path)

        # Write report
        write_threshold_report(thresholds, args.outdir / "threshold_report.txt")


if __name__ == "__main__":
    main()
