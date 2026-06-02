#!/usr/bin/env python3
"""
author: Yasas Wijesekara (yasas.wijesekara@uni-greifswald.de)

filters and formats the pileup output
"""
import polars as pl
import click
from taxopy.core import TaxDb
from taxopy import Taxon
from skadi.color_logger import logger


def get_acc2taxid_map(i:str) -> dict:
    """load accession2taxid and return a map"""
    acc2taxid = pl.read_csv(i,
            separator="\t",
            columns=["accession", "taxid"])
    
    return dict(zip(acc2taxid["accession"], acc2taxid["taxid"])) 

@click.command(
        short_help="summarize results of read annotation workflow"
)
@click.option(
    "-i", "--input",
    type=click.Path(exists=True, resolve_path=True),
    required=True,
    help="Input file"
)
@click.option(
    "-d", "--dbdir",
    type=click.Path(exists=True, resolve_path=True),
    required=True,
    help="skadi database root"
)
@click.option(
    "-o", "--output",
    type=click.Path(resolve_path=True),
    required=True,
    help="Output file"
)
@click.option(
    "-cp",
    type=float,
    required=False,
    default=60.0,
    help="Percentage of convered bases"
)
@click.option(
    "-af",
    type=float,
    required=False,
    default=1,
    help="Average fold coverage"
)
@click.option(
    "-mtr",
    type=float,
    required=False,
    default=100,
    help="Minimum total reads"
)
def main(input, dbdir, cp, af, mtr, output):
    """Filter and summarize pileup.sh output."""
    logger.info("summarizing pileup output")
    acc2taxid = get_acc2taxid_map(f"{dbdir}/VMR_latest/virus_genome.accession2taxid")

    taxdb = TaxDb(
        nodes_dmp=f"{dbdir}/ictv-taxdump/nodes.dmp",
        names_dmp=f"{dbdir}/ictv-taxdump/names.dmp",
        merged_dmp=f"{dbdir}/ictv-taxdump/merged.dmp",
    )

    def taxon_str(taxid: int) -> str:
        """Return a string representation of the Taxon lineage for a taxid."""
        if taxid is None:
            return ""
        try:
            return str(Taxon(taxid=int(taxid), taxdb=taxdb))
        except Exception:
            return ""

    df = pl.read_csv(input, separator="\t")

    df = (
        df.with_columns(
            total_reads=pl.col("Plus_reads") + pl.col("Minus_reads"),
            accession=pl.col("#ID").str.extract(r"^([A-Z]+\d+)", 1),
        )
        .with_columns(
            taxid=pl.col("accession").replace(acc2taxid),
        )
        .with_columns(
            lineage=pl.col("taxid").map_elements(taxon_str, return_dtype=pl.Utf8),
        )
        .filter(
            (pl.col("Covered_percent") >= cp)
            & (pl.col("Avg_fold") >= af)
            & (pl.col("total_reads") >= mtr)
        )
        .sort("Avg_fold", descending=True)
    )

    # Normalize column names last (so earlier references keep working)
    df = df.rename({c: c.lower() for c in df.columns})
    df.write_csv(output, separator="\t")

if __name__ == "__main__":
    main()
