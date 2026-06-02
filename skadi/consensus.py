"""Multi-method consensus taxonomy assignment.

Replaces the all-or-nothing ANI → AAI → API cascade with a weighted
voting scheme that combines evidence from all three methods.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple
import numpy as np
import polars as pl


# Default method weights derived from empirical accuracy benchmarks.
# ANI is most reliable for species/genus, AAI for family+, API as fallback.
DEFAULT_WEIGHTS = {
    "ani": {"species": 0.70, "genus": 0.65, "family": 0.40, "order": 0.30, "class": 0.25, "phylum": 0.20, "kingdom": 0.15},
    "aai": {"species": 0.20, "genus": 0.30, "family": 0.45, "order": 0.40, "class": 0.35, "phylum": 0.30, "kingdom": 0.25},
    "api": {"species": 0.10, "genus": 0.05, "family": 0.15, "order": 0.30, "class": 0.40, "phylum": 0.50, "kingdom": 0.60},
}

# Rank order for traversal (fine to coarse)
RANKS = ["species", "genus", "family", "order", "class", "phylum", "kingdom"]


def _score_to_confidence(score: float, threshold: float, steepness: float = 10.0) -> float:
    """Convert a raw score to a confidence using a sigmoid-like function.

    Args:
        score: Raw score (e.g., tani, taai, tapi).
        threshold: The rank-specific threshold.
        steepness: How sharply confidence rises above the threshold.

    Returns:
        Confidence in [0, 1].
    """
    # Logistic scaling: confidence = 1 / (1 + exp(-steepness * (score - threshold)))
    # At score == threshold, confidence = 0.5
    return float(1.0 / (1.0 + np.exp(-steepness * (score - threshold))))


def _get_rank_taxid_from_lineage(taxid: int, rank: str, taxdb) -> Optional[int]:
    """Get the taxid for a specific rank from a lineage."""
    try:
        from taxopy import Taxon
        taxon = Taxon(taxid, taxdb=taxdb)
        return taxon.rank_taxid_dictionary.get(rank)
    except Exception:
        return None


def build_consensus_assignment(
    ani_df: Optional[pl.DataFrame],
    aai_df: Optional[pl.DataFrame],
    api_df: Optional[pl.DataFrame],
    thresholds: Dict[str, float],
    taxdb,
    weights: Optional[Dict[str, Dict[str, float]]] = None,
    min_confidence: float = 0.5,
) -> pl.DataFrame:
    """Build consensus taxonomy assignment from ANI, AAI, and API results.

    Args:
        ani_df: DataFrame with columns [query, taxid, tani, ani, qcov, scov, ...].
        aai_df: DataFrame with columns [seqid, taxid, taai, aai, qcov, ...].
        api_df: DataFrame with columns [seqid, taxid, tapi, api, qcov, ...].
        thresholds: Dict of rank → threshold (e.g., {"species": 0.81, "genus": 0.49}).
        taxdb: Taxopy TaxDb object.
        weights: Optional custom method weights per rank.
        min_confidence: Minimum confidence required for assignment.

    Returns:
        DataFrame with columns [SequenceID, rank_taxid, rank_name, confidence, method].
    """
    weights = weights or DEFAULT_WEIGHTS

    # Collect all sequence IDs
    all_ids: set[str] = set()
    for df in (ani_df, aai_df, api_df):
        if df is not None and not df.is_empty():
            id_col = "query" if "query" in df.columns else "seqid"
            all_ids.update(df[id_col].to_list())

    rows = []
    for seqid in all_ids:
        # Gather votes from all methods
        votes: Dict[str, List[Tuple[int, float, str]]] = {rank: [] for rank in RANKS}

        # ANI votes
        if ani_df is not None and not ani_df.is_empty() and "query" in ani_df.columns:
            ani_row = ani_df.filter(pl.col("query") == seqid)
            if not ani_row.is_empty():
                taxid = ani_row["taxid"][0]
                tani = ani_row["tani"][0] if "tani" in ani_row.columns else ani_row["ani"][0]
                for rank in RANKS:
                    rank_taxid = _get_rank_taxid_from_lineage(taxid, rank, taxdb)
                    if rank_taxid is not None:
                        thresh = thresholds.get(f"tani{rank[0]}", thresholds.get(rank, 0.3))
                        conf = _score_to_confidence(tani, thresh)
                        votes[rank].append((rank_taxid, conf, "ani"))

        # AAI votes
        if aai_df is not None and not aai_df.is_empty() and "seqid" in aai_df.columns:
            aai_row = aai_df.filter(pl.col("seqid") == seqid)
            if not aai_row.is_empty():
                taxid = aai_row["taxid"][0]
                taai = aai_row["taai"][0] if "taai" in aai_row.columns else aai_row["aai"][0]
                for rank in RANKS:
                    rank_taxid = _get_rank_taxid_from_lineage(taxid, rank, taxdb)
                    if rank_taxid is not None:
                        thresh = thresholds.get(f"taai{rank[0]}", thresholds.get(rank, 0.3))
                        conf = _score_to_confidence(taai, thresh)
                        votes[rank].append((rank_taxid, conf, "aai"))

        # API votes
        if api_df is not None and not api_df.is_empty() and "seqid" in api_df.columns:
            api_row = api_df.filter(pl.col("seqid") == seqid)
            if not api_row.is_empty():
                taxid = api_row["taxid"][0]
                tapi = api_row["tapi"][0] if "tapi" in api_row.columns else api_row["api"][0]
                for rank in RANKS:
                    rank_taxid = _get_rank_taxid_from_lineage(taxid, rank, taxdb)
                    if rank_taxid is not None:
                        thresh = thresholds.get(f"tapi{rank[0]}", thresholds.get(rank, 0.15))
                        conf = _score_to_confidence(tapi, thresh)
                        votes[rank].append((rank_taxid, conf, "api"))

        # Weighted voting per rank
        best_rank = None
        best_taxid = None
        best_conf = 0.0
        best_methods = []

        for rank in RANKS:
            if not votes[rank]:
                continue

            # Weight confidence by method reliability at this rank
            weighted_votes: Dict[int, List[Tuple[float, str]]] = {}
            for taxid, conf, method in votes[rank]:
                w = weights.get(method, {}).get(rank, 0.33)
                weighted_conf = conf * w
                if taxid not in weighted_votes:
                    weighted_votes[taxid] = []
                weighted_votes[taxid].append((weighted_conf, method))

            # Sum weighted confidences per taxid
            for taxid, conf_methods in weighted_votes.items():
                total_conf = sum(cm[0] for cm in conf_methods)
                if total_conf > best_conf:
                    best_conf = total_conf
                    best_taxid = taxid
                    best_rank = rank
                    best_methods = [cm[1] for cm in conf_methods]

        # Only assign if confidence exceeds minimum
        if best_rank is not None and best_conf >= min_confidence:
            rows.append({
                "SequenceID": seqid,
                "rank": best_rank,
                "rank_taxid": best_taxid,
                "confidence": round(best_conf, 4),
                "methods": ",".join(sorted(set(best_methods))),
            })
        else:
            # Unassigned — could flag for manual review
            rows.append({
                "SequenceID": seqid,
                "rank": "unassigned",
                "rank_taxid": None,
                "confidence": round(best_conf, 4) if best_conf > 0 else 0.0,
                "methods": "",
            })

    return pl.DataFrame(rows) if rows else pl.DataFrame(
        {"SequenceID": [], "rank": [], "rank_taxid": [], "confidence": [], "methods": []}
    )
