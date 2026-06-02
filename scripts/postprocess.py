#!/usr/bin/env python
"""
author: Yasas Wijesekara (yasas.wijesekara@uni-greifswald.de)

postprocesses genome.m8, prot_dfein.m8 and prof_dfile.m8 files by calculating
ani, aai, and api and summarizes the summarizes the taxonomy predictions
to a single .tsv file

"""

import polars as pl
import click
from pathlib import Path
from taxopy.core import TaxDb
from taxopy import Taxon
from skadi.color_logger import logger


def trim_lineage(taxid: int, taxdb: TaxDb, taxomic_level: str = "species") -> int:
    # trims lineages to the level
    taxon = Taxon(taxid, taxdb)
    trimmed = taxon.rank_taxid_dictionary.get(taxomic_level)
    if trimmed is None:
        return taxid
    return trimmed


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


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.argument("db_dir", type=click.Path(exists=True, file_okay=False, path_type=str))
@click.argument("nuc", type=click.Path(exists=True, dir_okay=False, path_type=str))
@click.argument("prot", type=click.Path(exists=True, dir_okay=False, path_type=str))
@click.argument("prof", type=click.Path(exists=True, dir_okay=False, path_type=str))
@click.argument("outfile", type=click.Path(dir_okay=False, path_type=str))
# These were argv[6] and argv[7] but “should be moved to another place” → make them options:
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
def main(db_dir: str, nuc: str, prot: str, prof: str, outfile: str, **kwargs):
    # fasta = pyfastx.Fasta(FASTA)
    # headers = list(fasta.keys())
    # load the tax db
    taxdb = TaxDb(
        nodes_dmp=f"{db_dir}/ictv-taxdump/nodes.dmp",
        names_dmp=f"{db_dir}/ictv-taxdump/names.dmp",
        merged_dmp=f"{db_dir}/ictv-taxdump/merged.dmp",
    )

    def get_taxon(taxid: int) -> Taxon:
        """Instantiate a Taxon using the global `taxdb`."""
        return Taxon(taxid=taxid, taxdb=taxdb)

    def get_rank_taxid(taxid: int, rank: str) -> int:
        """Return the taxid that corresponds to *rank* (e.g. 'species', 'genus')."""
        return get_taxon(taxid).rank_taxid_dictionary.get(rank, 0)

    dfs = []
    matched = []
    try:
        # ani based filtering
        nuc_df = pl.read_csv(nuc, separator="\t")

        nuc_df = nuc_df.with_columns(pl.lit("ani").alias("Method")).rename(
            {"ani": "Score", "qlen": "Seqlen", "query": "SequenceID"}
        )

        df_species = nuc_df.filter(pl.col("tani") >= kwargs["tanis"]).with_columns(
            # taxid of the species that this row belongs to
            pl.col("taxid")
            .map_elements(lambda x: get_rank_taxid(x, "species"))
            .alias("rank_taxid"),
            # full lineage (Taxon object)
            pl.col("taxid")
            .map_elements(lambda x: str(get_taxon(x)))
            .alias("taxlineage"),
            # level label
            pl.lit("species").alias("level"),
        )
        df_genus = nuc_df.filter(
            (pl.col("tani") >= kwargs["tanig"]) & (pl.col("tani") < kwargs["tanis"])
        ).with_columns(
            # taxid of the genus that this row belongs to
            pl.col("taxid")
            .map_elements(lambda x: get_rank_taxid(x, "genus"))
            .alias("rank_taxid"),
            # full lineage (Taxon object)
            pl.col("taxid")
            .map_elements(lambda x: str(get_taxon(x)))
            .alias("taxlineage"),
            # level label
            pl.lit("genus").alias("level"),
        )

        # ------------------------------------------------------------------
        # Concatenate the two result DataFrames
        # ------------------------------------------------------------------
        nuc_df = pl.concat([df_species, df_genus])
        matched = nuc_df["SequenceID"].to_list()
        if not nuc_df.is_empty():
            dfs.append(nuc_df)
    except Exception as e:
        logger.info(f"nucleotide level results were not added because {e}")

    try:
        # aai based filtering
        prot_df = pl.read_csv(prot, separator="\t")

        prot_df = prot_df.with_columns(pl.lit("aai").alias("Method")).rename(
            {"taai": "Score", "qseqlen": "Seqlen", "seqid": "SequenceID"}
        )
        prot_df = prot_df.filter(~pl.col("SequenceID").is_in(matched))
        matched.extend(prot_df["SequenceID"].to_list())
        if not prot_df.is_empty():
            dfs.append(prot_df)

    except Exception:
        logger.info("no prot level results to merge")

    try:
        # api based filtering
        prof_df = pl.read_csv(prof, separator="\t")
        prof_df = prof_df.with_columns(pl.lit("api").alias("Method")).rename(
            {"tapi": "Score", "qseqlen": "Seqlen", "seqid": "SequenceID"}
        )

        prof_df = prof_df.filter(~pl.col("SequenceID").is_in(matched))
        matched.extend(prof_df["SequenceID"].to_list())
        if not prof_df.is_empty():
            dfs.append(prof_df)

    except Exception:
        logger.info("no profile level results to merge")

    # write a ictv taxonomy challange formatted file - this will be removed later

    if len(dfs) > 0:
        # logger.info(nuc_df.columns, keys)
        pl.concat([i.with_columns(*n).select(keys) for i in dfs]).write_csv(
            outfile.rstrip(".tsv") + "_ictv.csv", separator=","
        )

        # write a file with more information
        pl.concat([i.with_columns(*n).select(keys_full) for i in dfs]).write_csv(
            outfile, separator="\t"
        )
    else:
        Path(outfile).touch()


if __name__ == "__main__":
    main()

