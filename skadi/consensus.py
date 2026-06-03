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

# Threshold key prefixes used in config.yaml
_THRESH_PREFIXES = {"ani": "tani", "aai": "taai", "api": "tapi"}


def _resolve_threshold(thresholds: Dict[str, float], method: str, rank: str) -> float:
    """Resolve threshold for a method+rank from the threshold dict.

    Tries multiple key formats:
      - "tanis" (legacy prefix+rank_initial)
      - "tani_species" (prefix_rank)
      - "species" (rank name only)
      - fallback default
    """
    prefix = _THRESH_PREFIXES.get(method, method)
    # Try legacy format: tanis, tanig, etc.
    legacy_key = f"{prefix}{rank[0]}"
    if legacy_key in thresholds:
        return thresholds[legacy_key]
    # Try prefix_rank format
    modern_key = f"{prefix}_{rank}"
    if modern_key in thresholds:
        return thresholds[modern_key]
    # Try rank-only
    if rank in thresholds:
        return thresholds[rank]
    # Fallback defaults by method
    defaults = {"ani": 0.3, "aai": 0.3, "api": 0.15}
    return defaults.get(method, 0.3)


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


def _get_best_hit(df: pl.DataFrame, id_col: str, score_col: str) -> Optional[pl.DataFrame]:
    """Get the best hit (highest score) for each sequence ID."""
    if df is None or df.is_empty():
        return None
    return df.sort(score_col, descending=True).group_by(id_col).first()


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
        DataFrame with columns [SequenceID, rank, rank_taxid, confidence, methods,
                                taxlineage, Score, level]. Compatible with postprocess.py.
    """
    weights = weights or DEFAULT_WEIGHTS

    # Get best hit per sequence for each method
    ani_best = _get_best_hit(ani_df, "query", "tani") if ani_df is not None else None
    aai_best = _get_best_hit(aai_df, "seqid", "taai") if aai_df is not None else None
    api_best = _get_best_hit(api_df, "seqid", "tapi") if api_df is not None else None

    # Collect all sequence IDs
    all_ids: set[str] = set()
    for df in (ani_best, aai_best, api_best):
        if df is not None and not df.is_empty():
            id_col = "query" if "query" in df.columns else "seqid"
            all_ids.update(df[id_col].to_list())

    rows = []
    for seqid in all_ids:
        # Gather votes from all methods
        votes: Dict[str, List[Tuple[int, float, str, float]]] = {rank: [] for rank in RANKS}

        # ANI votes (best hit)
        if ani_best is not None and not ani_best.is_empty() and "query" in ani_best.columns:
            ani_row = ani_best.filter(pl.col("query") == seqid)
            if not ani_row.is_empty():
                taxid = int(ani_row["taxid"][0])
                tani = float(ani_row["tani"][0]) if "tani" in ani_row.columns else float(ani_row["ani"][0])
                for rank in RANKS:
                    rank_taxid = _get_rank_taxid_from_lineage(taxid, rank, taxdb)
                    if rank_taxid is not None:
                        thresh = _resolve_threshold(thresholds, "ani", rank)
                        conf = _score_to_confidence(tani, thresh)
                        votes[rank].append((rank_taxid, conf, "ani", tani))

        # AAI votes (best hit)
        if aai_best is not None and not aai_best.is_empty() and "seqid" in aai_best.columns:
            aai_row = aai_best.filter(pl.col("seqid") == seqid)
            if not aai_row.is_empty():
                taxid = int(aai_row["taxid"][0])
                taai = float(aai_row["taai"][0]) if "taai" in aai_best.columns else float(aai_row["aai"][0])
                for rank in RANKS:
                    rank_taxid = _get_rank_taxid_from_lineage(taxid, rank, taxdb)
                    if rank_taxid is not None:
                        thresh = _resolve_threshold(thresholds, "aai", rank)
                        conf = _score_to_confidence(taai, thresh)
                        votes[rank].append((rank_taxid, conf, "aai", taai))

        # API votes (best hit)
        if api_best is not None and not api_best.is_empty() and "seqid" in api_best.columns:
            api_row = api_best.filter(pl.col("seqid") == seqid)
            if not api_row.is_empty():
                taxid = int(api_row["taxid"][0])
                tapi = float(api_row["tapi"][0]) if "tapi" in api_best.columns else float(api_row["api"][0])
                for rank in RANKS:
                    rank_taxid = _get_rank_taxid_from_lineage(taxid, rank, taxdb)
                    if rank_taxid is not None:
                        thresh = _resolve_threshold(thresholds, "api", rank)
                        conf = _score_to_confidence(tapi, thresh)
                        votes[rank].append((rank_taxid, conf, "api", tapi))

        # Weighted voting per rank
        best_rank = None
        best_taxid = None
        best_conf = 0.0
        best_methods = []
        best_score = 0.0

        for rank in RANKS:
            if not votes[rank]:
                continue

            # Weight confidence by method reliability at this rank
            weighted_votes: Dict[int, List[Tuple[float, str, float]]] = {}
            for taxid, conf, method, raw_score in votes[rank]:
                w = weights.get(method, {}).get(rank, 0.33)
                weighted_conf = conf * w
                if taxid not in weighted_votes:
                    weighted_votes[taxid] = []
                weighted_votes[taxid].append((weighted_conf, method, raw_score))

            # Sum weighted confidences per taxid
            for taxid, conf_methods in weighted_votes.items():
                total_conf = sum(cm[0] for cm in conf_methods)
                if total_conf > best_conf:
                    best_conf = total_conf
                    best_taxid = taxid
                    best_rank = rank
                    best_methods = [cm[1] for cm in conf_methods]
                    best_score = max(cm[2] for cm in conf_methods)

        # Build lineage string for the best taxid
        taxlineage = ""
        if best_taxid is not None:
            try:
                from taxopy import Taxon
                taxlineage = str(Taxon(best_taxid, taxdb=taxdb))
            except Exception:
                taxlineage = ""

        # Only assign if confidence exceeds minimum
        if best_rank is not None and best_conf >= min_confidence:
            rows.append({
                "SequenceID": seqid,
                "rank": best_rank,
                "rank_taxid": best_taxid,
                "confidence": round(best_conf, 4),
                "methods": ",".join(sorted(set(best_methods))),
                "taxlineage": taxlineage,
                "Score": round(best_score, 4),
                "level": best_rank,
            })
        else:
            # Unassigned — could flag for manual review
            rows.append({
                "SequenceID": seqid,
                "rank": "unassigned",
                "rank_taxid": None,
                "confidence": round(best_conf, 4) if best_conf > 0 else 0.0,
                "methods": "",
                "taxlineage": "",
                "Score": 0.0,
                "level": "unassigned",
            })

    return pl.DataFrame(rows) if rows else pl.DataFrame(
        {"SequenceID": [], "rank": [], "rank_taxid": [], "confidence": [],
         "methods": [], "taxlineage": [], "Score": [], "level": []}
    )
