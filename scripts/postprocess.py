#!/usr/bin/env python
"""
author: Yasas Wijesekara (yasas.wijesekara@uni-greifswald.de)

postprocesses genome.m8, prot_dfein.m8 and prof_dfile.m8 files by calculating
ani, aai, and api and summarizes the summarizes the taxonomy predictions
to a single .tsv file

Supports two assignment methods:
  - cascade (default): ANI → AAI → API tiered fallback
  - consensus: weighted voting across all three methods

"""

from typing import Dict, List, Optional
import json
import polars as pl
import click
from pathlib import Path
from taxopy.core import TaxDb
from taxopy import Taxon
from skadi.color_logger import logger
from skadi.consensus import build_consensus_assignment


keys = [
    "SequenceID",
    "Realm (-viria)",
    "Realm_score",
    "Subrealm (-vira)",
    "Subrealm_score",
    "Kingdom (-virae)",
    "Kingdom_score",
    "Subkingdom (-virites)",
    "Subkingdom_score",
    "Phylum (-viricota)",
    "Phylum_score",
    "Subphylum (-viricotina)",
    "Subphylum_score",
    "Class (-viricetes)",
    "Class_score",
    "Subclass (-viricetidae)",
    "Subclass_score",
    "Order (-virales)",
    "Order_score",
    "Suborder (-virineae)",
    "Suborder_score",
    "Family (-viridae)",
    "Family_score",
    "Subfamily (-virinae)",
    "Subfamily_score",
    "Genus (-virus)",
    "Genus_score",
    "Subgenus (-virus)",
    "Subgenus_score",
    "Species (binomial)",
    "Species_score",
]

keys_full = [
    "SequenceID",
    "Seqlen",
    "Score",
    "Method",
    "Realm (-viria)",
    "Realm_score",
    "Subrealm (-vira)",
    "Subrealm_score",
    "Kingdom (-virae)",
    "Kingdom_score",
    "Subkingdom (-virites)",
    "Subkingdom_score",
    "Phylum (-viricota)",
    "Phylum_score",
    "Subphylum (-viricotina)",
    "Subphylum_score",
    "Class (-viricetes)",
    "Class_score",
    "Subclass (-viricetidae)",
    "Subclass_score",
    "Order (-virales)",
    "Order_score",
    "Suborder (-virineae)",
    "Suborder_score",
    "Family (-viridae)",
    "Family_score",
    "Subfamily (-virinae)",
    "Subfamily_score",
    "Genus (-virus)",
    "Genus_score",
    "Subgenus (-virus)",
    "Subgenus_score",
    "Species (binomial)",
    "Species_score",
]
n = [
    pl.col("taxlineage").str.extract(r"_([A-Za-z]+viria);?", 1).alias("Realm (-viria)"),
    pl.lit(None).alias("Realm_score"),
    pl.col("taxlineage").str.extract(r"_([A-Za-z]+vira);", 1).alias("Subrealm (-vira)"),
    pl.lit(None).alias("Subrealm_score"),
    pl.col("taxlineage")
    .str.extract(r"_([A-Za-z]+virae);?", 1)
    .alias("Kingdom (-virae)"),
    pl.lit(None).alias("Kingdom_score"),
    pl.col("taxlineage")
    .str.extract(r"_([A-Za-z]+virites);?", 1)
    .alias("Subkingdom (-virites)"),
    pl.lit(None).alias("Subkingdom_score"),
    pl.col("taxlineage")
    .str.extract(r"_([A-Za-z]+viricota);?", 1)
    .alias("Phylum (-viricota)"),
    pl.lit(None).alias("Phylum_score"),
    pl.col("taxlineage")
    .str.extract(r"_([A-Za-z]+viricotina);?", 1)
    .alias("Subphylum (-viricotina)"),
    pl.lit(None).alias("Subphylum_score"),
    pl.col("taxlineage")
    .str.extract(r"_([A-Za-z]+viricetes);?", 1)
    .alias("Class (-viricetes)"),
    pl.lit(None).alias("Class_score"),
    pl.col("taxlineage")
    .str.extract(r"_([A-Za-z]+viricetidae);?", 1)
    .alias("Subclass (-viricetidae)"),
    pl.lit(None).alias("Subclass_score"),
    pl.col("taxlineage")
    .str.extract(r"_([A-Za-z]+virales);?", 1)
    .alias("Order (-virales)"),
    pl.lit(None).alias("Order_score"),
    pl.col("taxlineage")
    .str.extract(r"_([A-Za-z]+virineae);?", 1)
    .alias("Suborder (-virineae)"),
    pl.lit(None).alias("Suborder_score"),
    pl.col("taxlineage")
    .str.extract(r"_([A-Za-z]+viridae);?", 1)
    .alias("Family (-viridae)"),
    pl.lit(None).alias("Family_score"),
    pl.col("taxlineage")
    .str.extract(r"_([A-Za-z]+virinae);?", 1)
    .alias("Subfamily (-virinae)"),
    pl.lit(None).alias("Subfamily_score"),
    pl.col("taxlineage").str.extract(r"_([A-Za-z]+virus);?", 1).alias("Genus (-virus)"),
    pl.lit(None).alias("Genus_score"),
    pl.col("taxlineage")
    .str.extract(r"_[A-Za-z]+virus;-_([A-Za-z]+virus);?", 1)
    .alias("Subgenus (-virus)"),
    pl.lit(None).alias("Subgenus_score"),
    pl.col("taxlineage")
    .str.extract(r"_([A-Za-z]+(?:\s[A-Za-z0-9-]+)+);?", 1)
    .alias("Species (binomial)"),
    pl.lit(None).alias("Species_score"),
]


def _build_taxon_lookup(taxdb: TaxDb, taxids: List[int]) -> Dict[int, Taxon]:
    """Pre-build Taxon objects for a list of taxids."""
    return {tid: Taxon(tid, taxdb=taxdb) for tid in set(taxids) if tid}


def _build_rank_lookup(taxdb: TaxDb, taxids: List[int], rank: str) -> Dict[int, int]:
    """Pre-build taxid → rank_taxid mapping for vectorized operations."""
    lookup = {}
    for tid in set(taxids):
        if tid:
            try:
                lookup[tid] = Taxon(tid, taxdb=taxdb).rank_taxid_dictionary.get(rank, 0)
            except Exception:
                lookup[tid] = 0
    return lookup


def _build_lineage_lookup(taxdb: TaxDb, taxids: List[int]) -> Dict[int, str]:
    """Pre-build taxid → lineage string mapping for vectorized operations."""
    lookup = {}
    for tid in set(taxids):
        if tid:
            try:
                lookup[tid] = str(Taxon(tid, taxdb=taxdb))
            except Exception:
                lookup[tid] = ""
    return lookup


def _apply_taxon_columns(df: pl.DataFrame, taxdb: TaxDb, rank: str) -> pl.DataFrame:
    """Vectorized addition of rank_taxid and taxlineage columns."""
    unique_taxids = [t for t in df["taxid"].unique().to_list() if t is not None]
    rank_map = _build_rank_lookup(taxdb, unique_taxids, rank)
    lineage_map = _build_lineage_lookup(taxdb, unique_taxids)
    return df.with_columns(
        pl.col("taxid").replace_strict(rank_map, default=0).alias("rank_taxid"),
        pl.col("taxid").replace_strict(lineage_map, default="").alias("taxlineage"),
        pl.lit(rank).alias("level"),
    )


def _load_thresholds(thresholds_file: str) -> Dict:
    """Load per-group thresholds from JSON file.

    Supports both old flat format and new nested format (ani/aai/api).

    Returns:
        Dict with structure:
        {
            "ani": {"global": {...}, "species": {"family": {...}}, "genus": {"family": {...}}},
            "aai": {"global": {...}, "family": {"order": {...}}, ...},
            "api": {"global": {...}, "family": {"order": {...}}, ...},
        }
        Or falls back to old flat format.
    """
    with open(thresholds_file) as f:
        data = json.load(f)

    # Detect format: new format has "ani" or "aai" or "api" as top-level keys
    if any(k in data for k in ["ani", "aai", "api"]):
        return data

    # Old flat format: wrap as "ani" for backward compatibility
    return {"ani": data}


def _run_cascade(
    db_dir: str, nuc: str, prot: str, prof: str, taxdb: TaxDb, kwargs: dict
) -> List[pl.DataFrame]:
    """Original ANI → AAI → API cascade logic.

    Supports per-group thresholds via --thresholds-file.
    """
    dfs = []
    matched = []

    # Load thresholds if provided
    thresholds_data = None
    ani_thresholds = None
    aai_thresholds = None
    api_thresholds = None
    thresholds_file = kwargs.get("thresholds_file")
    if thresholds_file and Path(thresholds_file).exists():
        thresholds_data = _load_thresholds(thresholds_file)
        ani_thresholds = thresholds_data.get("ani", {})
        aai_thresholds = thresholds_data.get("aai", {})
        api_thresholds = thresholds_data.get("api", {})
        logger.info("Loaded thresholds from %s", thresholds_file)

    try:
        nuc_df = pl.read_csv(nuc, separator="\t")
        nuc_df = nuc_df.with_columns(pl.lit("ani").alias("Method")).rename(
            {"ani": "Score", "qlen": "Seqlen", "query": "SequenceID"}
        )

        # Get ANI thresholds (species + genus)
        ani_global = ani_thresholds.get("global", {}) if ani_thresholds else {}
        tanis = ani_global.get("species", {}).get("threshold", kwargs["tanis"])
        tanig = ani_global.get("genus", {}).get("threshold", kwargs["tanig"])

        # Add family_taxid for per-family threshold lookup
        unique_taxids = [t for t in nuc_df["taxid"].unique().to_list() if t is not None]
        family_map = _build_rank_lookup(taxdb, unique_taxids, "family")
        nuc_df = nuc_df.with_columns(
            pl.col("taxid").replace_strict(family_map, default=0).alias("family_taxid")
        )

        # Build per-family ANI thresholds dataframe
        if ani_thresholds:
            family_taxids = set()
            for rank in ["species", "genus"]:
                if rank in ani_thresholds and "family" in ani_thresholds[rank]:
                    family_taxids.update(ani_thresholds[rank]["family"].keys())

            rows = []
            for taxid in family_taxids:
                sp = ani_thresholds.get("species", {}).get("family", {}).get(taxid, {}).get("threshold")
                gn = ani_thresholds.get("genus", {}).get("family", {}).get(taxid, {}).get("threshold")
                rows.append({"family_taxid": int(taxid), "species_thresh": sp, "genus_thresh": gn})

            if rows:
                family_thresh_df = pl.DataFrame(rows)
                nuc_df = nuc_df.join(family_thresh_df, on="family_taxid", how="left")
                nuc_df = nuc_df.with_columns(
                    pl.col("species_thresh").fill_null(tanis),
                    pl.col("genus_thresh").fill_null(tanig),
                )
            else:
                nuc_df = nuc_df.with_columns(
                    pl.lit(tanis).alias("species_thresh"),
                    pl.lit(tanig).alias("genus_thresh"),
                )
        else:
            nuc_df = nuc_df.with_columns(
                pl.lit(tanis).alias("species_thresh"),
                pl.lit(tanig).alias("genus_thresh"),
            )

        # Filter using per-family thresholds
        df_species = nuc_df.filter(pl.col("tani") >= pl.col("species_thresh"))
        df_species = _apply_taxon_columns(df_species, taxdb, "species")

        df_genus = nuc_df.filter(
            (pl.col("tani") >= pl.col("genus_thresh")) & (pl.col("tani") < pl.col("species_thresh"))
        )
        df_genus = _apply_taxon_columns(df_genus, taxdb, "genus")

        nuc_df = pl.concat([df_species, df_genus])
        matched = nuc_df["SequenceID"].to_list()
        if not nuc_df.is_empty():
            dfs.append(nuc_df)
    except (pl.exceptions.ComputeError, FileNotFoundError, ValueError) as e:
        logger.info(f"nucleotide level results were not added because {e}")

    try:
        prot_df = pl.read_csv(prot, separator="\t")
        prot_df = prot_df.with_columns(pl.lit("aai").alias("Method")).rename(
            {"taai": "Score", "qseqlen": "Seqlen", "seqid": "SequenceID"}
        )
        prot_df = prot_df.filter(~pl.col("SequenceID").is_in(matched))
        matched.extend(prot_df["SequenceID"].to_list())
        if not prot_df.is_empty():
            dfs.append(prot_df)
    except (pl.exceptions.ComputeError, FileNotFoundError, ValueError):
        logger.info("no prot level results to merge")

    try:
        prof_df = pl.read_csv(prof, separator="\t")
        prof_df = prof_df.with_columns(pl.lit("api").alias("Method")).rename(
            {"tapi": "Score", "qseqlen": "Seqlen", "seqid": "SequenceID"}
        )
        prof_df = prof_df.filter(~pl.col("SequenceID").is_in(matched))
        matched.extend(prof_df["SequenceID"].to_list())
        if not prof_df.is_empty():
            dfs.append(prof_df)
    except (pl.exceptions.ComputeError, FileNotFoundError, ValueError):
        logger.info("no profile level results to merge")

    return dfs


def _run_consensus(
    db_dir: str, nuc: str, prot: str, prof: str, taxdb: TaxDb, kwargs: dict
) -> List[pl.DataFrame]:
    """Consensus-based taxonomy assignment using weighted voting."""
    ani_df = None
    aai_df = None
    api_df = None

    try:
        ani_df = pl.read_csv(nuc, separator="\t")
    except (pl.exceptions.ComputeError, FileNotFoundError, ValueError) as e:
        logger.info(f"nucleotide level results not available: {e}")

    try:
        aai_df = pl.read_csv(prot, separator="\t")
    except (pl.exceptions.ComputeError, FileNotFoundError, ValueError) as e:
        logger.info(f"protein level results not available: {e}")

    try:
        api_df = pl.read_csv(prof, separator="\t")
    except (pl.exceptions.ComputeError, FileNotFoundError, ValueError) as e:
        logger.info(f"profile level results not available: {e}")

    consensus = build_consensus_assignment(
        ani_df=ani_df,
        aai_df=aai_df,
        api_df=api_df,
        thresholds=kwargs,
        taxdb=taxdb,
        min_confidence=kwargs.get("min_confidence", 0.5),
    )

    if consensus.is_empty():
        return []

    # Convert consensus output to the format expected by downstream processing
    # Add Method column based on methods used
    consensus = consensus.with_columns(
        pl.col("methods").alias("Method"),
        pl.col("Score").alias("Score"),
        pl.col("SequenceID").alias("SequenceID"),
    )

    # For unassigned, we still want them in the output
    return [consensus]


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.argument("db_dir", type=click.Path(exists=True, file_okay=False, path_type=str))
@click.argument("nuc", type=click.Path(exists=True, dir_okay=False, path_type=str))
@click.argument("prot", type=click.Path(exists=True, dir_okay=False, path_type=str))
@click.argument("prof", type=click.Path(exists=True, dir_okay=False, path_type=str))
@click.argument("outfile", type=click.Path(dir_okay=False, path_type=str))
@click.option(
    "--method",
    type=click.Choice(["cascade", "consensus"], case_sensitive=False),
    default="cascade",
    help="Taxonomy assignment method: cascade (ANI→AAI→API tiered) or consensus (weighted voting).",
)
@click.option(
    "--min-confidence",
    type=float,
    default=0.5,
    help="Minimum confidence for consensus assignment (only used with --method consensus).",
)
@click.option(
    "--tapif",
    type=float,
    default=0.15,
    help="assign sequences above this tapi threshold to families",
    required=False,
)
@click.option(
    "--tapio",
    type=float,
    default=0.15,
    help="assign sequences above this tapi threshold to orders",
    required=False,
)
@click.option(
    "--tapic",
    type=float,
    default=0.15,
    help="assign sequences above this tapi threshold to classes",
    required=False,
)
@click.option(
    "--tapip",
    type=float,
    default=0.15,
    help="assign sequences above this tapi threshold to phyla",
    required=False,
)
@click.option(
    "--tapik",
    type=float,
    default=0.15,
    help="assign sequences above this tapi threshold to kingdoms",
    required=False,
)
@click.option(
    "--taaig",
    type=float,
    default=0.3,
    help="assign sequences above this taai threshold to genera",
    required=False,
)
@click.option(
    "--taaif",
    type=float,
    default=0.3,
    help="assign sequences above this taai threshold to families",
    required=False,
)
@click.option(
    "--taaio",
    type=float,
    default=0.3,
    help="assign sequences above this taai threshold to orders",
    required=False,
)
@click.option(
    "--taaic",
    type=float,
    default=0.3,
    help="assign sequences above this taai threshold to classes",
    required=False,
)
@click.option(
    "--taaip",
    type=float,
    default=0.3,
    help="assign sequences above this taai threshold to phyla",
    required=False,
)
@click.option(
    "--taaik",
    type=float,
    default=0.3,
    help="assign sequences above this taai threshold to kingdoms",
    required=False,
)
@click.option(
    "--tanis",
    type=float,
    default=0.81,
    help="assign sequences above this taai threshold to species",
    required=False,
)
@click.option(
    "--tanig",
    type=float,
    default=0.49,
    help="assign sequences above this taai threshold to genera",
    required=False,
)
@click.option(
    "--thresholds-file",
    type=click.Path(exists=True, dir_okay=False, path_type=str),
    default=None,
    help="JSON file with per-group thresholds from compute_pairwise_matrix.py",
)
def main(db_dir: str, nuc: str, prot: str, prof: str, outfile: str, method: str, **kwargs):
    taxdb = TaxDb(
        nodes_dmp=f"{db_dir}/ictv-taxdump/nodes.dmp",
        names_dmp=f"{db_dir}/ictv-taxdump/names.dmp",
        merged_dmp=f"{db_dir}/ictv-taxdump/merged.dmp",
    )

    if method == "consensus":
        dfs = _run_consensus(db_dir, nuc, prot, prof, taxdb, kwargs)
    else:
        dfs = _run_cascade(db_dir, nuc, prot, prof, taxdb, kwargs)

    if len(dfs) > 0:
        pl.concat([i.with_columns(*n).select(keys) for i in dfs]).write_csv(
            outfile.rstrip(".tsv") + "_ictv.csv", separator=","
        )
        pl.concat([i.with_columns(*n).select(keys_full) for i in dfs]).write_csv(
            outfile, separator="\t"
        )
    else:
        Path(outfile).touch()


if __name__ == "__main__":
    main()
