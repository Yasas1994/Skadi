#!/usr/bin/env python3
"""
author: Yasas Wijesekara (yasas.wijesekara@uni-greifswald.de)

This script adds lineage information to .m8 files with a taxid column.
Run after blast_genome rule.
"""


import logging
import re
from pathlib import Path

import click
import polars as pl
import taxopy

DEFAULT_IN_HEADER = (
    "query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,"
    "qstart,qend,tstart,tend,evalue,bits,taxid"
)

DEFAULT_OUT_HEADER = (
    "query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,"
    "qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage"
)

DEFAULT_TAXID = "taxid"


def setup_logging(level: str) -> logging.Logger:
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler()],
    )
    return logging.getLogger("skadi")


def search_pattern(x: str) -> str:
    m = re.search(r"(?<=\|)(\S+)(?=\|)", x)
    return m.group(1) if m else ""


@click.command(context_settings=dict(show_default=True))
@click.option(
    "-d",
    "--database-dir",
    type=click.Path(path_type=Path, exists=True, file_okay=False, dir_okay=True),
    required=True,
    help="Root database directory (contains ictv-taxdump/ and VMR_latest/).",
)
@click.option(
    "-i",
    "--input",
    "table_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="Input m8 file path.",
)
@click.option(
    "--in-header",
    default=DEFAULT_IN_HEADER,
    show_default=False,
    help="Comma-separated header for the input m8 file.",
)
@click.option(
    "--out-header",
    default=DEFAULT_OUT_HEADER,
    show_default=False,
    help="Comma-separated header for the output m8 file (kept for downstream compatibility).",
)
@click.option(
    "--taxid-col",
    default=DEFAULT_TAXID,
    show_default=True,
    help="Column name containing the taxid.",
)
@click.option(
    "--taxdump-dir",
    type=click.Path(path_type=Path, exists=True, file_okay=False, dir_okay=True),
    default=None,
    help="Override path to ictv-taxdump/ (defaults to <database-dir>/ictv-taxdump).",
)
@click.option(
    "--log-level",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"], case_sensitive=False),
    default="INFO",
    help="Logging verbosity.",
)
def main(
    database_dir: Path,
    table_path: Path,
    in_header: str,
    out_header: str,   # kept for compatibility; not required for writing since header=False
    taxid_col: str,
    taxdump_dir: Path | None,
    log_level: str,
) -> None:
    logger = setup_logging(log_level)

    database_dir = database_dir.resolve()
    taxdump_dir = (taxdump_dir or (database_dir / "ictv-taxdump")).resolve()
    table_path = table_path.resolve()
    table_path_out = table_path.with_suffix(table_path.suffix + ".tmp")

    # Validate required taxdump files
    required = [
        (taxdump_dir / "nodes.dmp", "nodes.dmp"),
        (taxdump_dir / "names.dmp", "names.dmp"),
        (taxdump_dir / "merged.dmp", "merged.dmp"),
    ]
    for p, desc in required:
        if not p.exists():
            raise click.ClickException(f"Required path for {desc} not found: {p}")

    logger.info("Using database dir: %s", database_dir)
    logger.info("Using taxdump dir:  %s", taxdump_dir)
    logger.info("Using m8 table:     %s", table_path)

    # Load taxonomy DB
    taxdb = taxopy.TaxDb(
        nodes_dmp=str(taxdump_dir / "nodes.dmp"),
        names_dmp=str(taxdump_dir / "names.dmp"),
        merged_dmp=str(taxdump_dir / "merged.dmp"),
    )
    logger.info("Loaded TaxDb with %d nodes.", len(taxdb.taxid2parent))

    # Read m8
    cols = [c.strip() for c in in_header.split(",") if c.strip()]
    table = pl.read_csv(
        table_path,
        separator="\t",
        has_header=False,
        new_columns=cols,
    )
    logger.info("Loaded m8: %s rows.", f"{table.height:,}")

    if taxid_col not in table.columns:
        raise click.ClickException(
            f"taxid-col '{taxid_col}' not found in columns: {', '.join(table.columns)}"
        )

    # Build Taxon cache only for unique non-null taxids (much faster + less memory)
    taxids = (
        table.select(pl.col(taxid_col))
        .unique()
        .drop_nulls()
        .to_series()
        .to_list()
    )

    taxid2taxon = {}
    for t in taxids:
        try:
            taxid2taxon[int(t)] = taxopy.Taxon(int(t), taxdb=taxdb)
        except Exception:
            # unknown/invalid taxid -> leave unmapped
            pass

    def get_name(x) -> str:
        tx = taxid2taxon.get(int(x)) if x is not None else None
        return tx.name if tx else "NA"

    def get_lineage(x) -> str:
        tx = taxid2taxon.get(int(x)) if x is not None else None
        return str(tx) if tx else "NA"

    # Transform + annotate
    table = table.with_columns(
        (pl.col("fident") / 100).alias("fident"),
        pl.col("target").map_elements(search_pattern, return_dtype=pl.String).alias("target"),
    ).with_columns(
        taxname=pl.col(taxid_col).map_elements(get_name, return_dtype=pl.String),
        taxlineage=pl.col(taxid_col).map_elements(get_lineage, return_dtype=pl.String),
    )

    table.write_csv(table_path_out, separator="\t", include_header=False)
    logger.info("Wrote: %s", table_path_out)


if __name__ == "__main__":
    main()
