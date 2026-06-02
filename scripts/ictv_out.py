#!/usr/bin/env python3
"""
author: Yasas Wijesekara (yasas.wijesekara@uni-greifswald.de)

Formats SKADI outputs to the ICTV taxonomy challenge result format.
"""

from pathlib import Path
import click
import pandas as pd
import taxopy

from skadi.color_logger import logger


@click.command(context_settings=dict(show_default=True))
@click.option(
    "-d",
    "--database-dir",
    type=click.Path(path_type=Path, exists=True, file_okay=False),
    required=True,
    help="Root database directory containing ictv-taxdump/.",
)
@click.option(
    "-i",
    "--input",
    "infile",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="SKADI result input file.",
)
@click.option(
    "-o",
    "--output",
    "outfile",
    type=click.Path(path_type=Path, dir_okay=False),
    required=True,
    help="Output file in ICTV taxonomy challenge format.",
)
def main(database_dir: Path, infile: Path, outfile: Path) -> None:
    database_dir = database_dir.resolve()
    infile = infile.resolve()
    outfile = outfile.resolve()

    # Load taxonomy database
    taxdb = taxopy.TaxDb(
        nodes_dmp=str(database_dir / "ictv-taxdump/nodes.dmp"),
        names_dmp=str(database_dir / "ictv-taxdump/names.dmp"),
        merged_dmp=str(database_dir / "ictv-taxdump/merged.dmp"),
    )

    logger.info("Loaded ICTV taxonomy database")

    results = pd.read_csv(infile)

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

    rows = []

    for _, row in results.iterrows():
        out_row = dict.fromkeys(keys)
        out_row["SequenceID"] = row["query"]

        taxon = taxopy.Taxon(row["taxid"], taxdb=taxdb)
        rank_map = taxon.rank_name_dictionary

        for rank, name in rank_map.items():
            for col in keys:
                if col.startswith(f"{rank.capitalize()} "):
                    out_row[col] = name

        rows.append(out_row)

    logger.info("Writing ICTV formatted output file: %s", outfile)
    pd.DataFrame(rows).to_csv(outfile, index=False)


if __name__ == "__main__":
    main()
