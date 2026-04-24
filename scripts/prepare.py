#!/usr/bin/env python3
"""
author: Yasas Wijesekara (yasas.wijesekara@uni-greifswald.de)

Prepares the genome and protein sequence files for MMseqs database creation.
This script needs refractoring!
"""

import re
from glob import glob
from pathlib import Path
from typing import List

import click
import numpy as np
import openpyxl
import pandas as pd
import taxopy
from Bio import GenBank, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from vcat.color_logger import logger


def extract_regions(string):
    return [list(map(int, i)) for i in re.findall(r"\((\d+)\.{1,3}(\d+)\)", str(string))]


def filter_genbank(ids):
    ids_ = [i for i in ids if "_" in i]
    return ids_ if ids_ else ids


def extract_from_hyperlink(hypstr):
    url_match = re.search(r'=HYPERLINK\("([^"]+)"', str(hypstr))
    url = url_match.group(1) if url_match else None
    text_match = re.search(r',"([^"]+)"\)', str(hypstr))
    display_text = text_match.group(1) if text_match else None
    return {"url": f"{url}", "display_text": f"{display_text}"}

def extract_species_name(text: str) -> str | None:
    """
    Extract species name from an Excel-style HYPERLINK formula like:
    =HYPERLINK("https://ictv.global/taxonomy/taxondetails?taxnode_id=202308643","Alphalipothrixvirus beppuense")

    If the hyperlink pattern is not found, return the input text stripped,
    assuming it already contains the species name.
    """
    if text is None:
        return None

    text = text.strip()

    pattern = r'^=HYPERLINK\(".*?","(.*?)"\)$'
    match = re.match(pattern, text)

    if match:
        return match.group(1).strip()

    return text if text else None


def get_taxid_from_species(species: str , taxdb: taxopy.TaxDb) -> List[int] | List[List[int]] | None:
    name = extract_species_name(species)
    if not name:
        return None
    taxids = taxopy.taxid_from_name(name, taxdb)
    return taxids


@click.command(context_settings=dict(show_default=True, help_option_names=["-h", "--help"]))
@click.option(
    "-d",
    "--database-dir",
    type=click.Path(path_type=Path, exists=True, file_okay=False),
    required=True,
    help="Database directory containing ictv-taxdump/.",
)
@click.option(
    "-x",
    "--xltable",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
    help="ICTV VMR Excel table (.xlsx).",
)
@click.option(
    "-g",
    "--gbdir",
    type=click.Path(path_type=Path, exists=True, file_okay=False),
    required=True,
    help="Directory containing downloaded GenBank (.gb) files.",
)
@click.option(
    "-s",
    "--seqdir",
    type=click.Path(path_type=Path, file_okay=False),
    required=True,
    help="Output directory for genomes.fna/proteins.faa and accession2taxid files.",
)
def main(database_dir: Path, xltable: Path, gbdir: Path, seqdir: Path) -> None:
    database_dir = database_dir.resolve()
    xltable = xltable.resolve()
    gbdir = gbdir.resolve()
    seqdir = seqdir.resolve()
    seqdir.mkdir(parents=True, exist_ok=True)

    taxdb = taxopy.TaxDb(
        nodes_dmp=str(database_dir / "ictv-taxdump/nodes.dmp"),
        names_dmp=str(database_dir / "ictv-taxdump/names.dmp"),
        merged_dmp=str(database_dir / "ictv-taxdump/merged.dmp"),
    )

    wb = openpyxl.load_workbook(xltable)
    ws = wb[wb.sheetnames[1]]
    data = ws.values
    columns = next(data)
    ictv = pd.DataFrame(data, columns=[str(i) for i in columns])

    ictv = ictv[ictv["Accessions Link"].notna()]

    ictv["Taxid"] = ictv["Species"].apply(lambda x: taxopy.taxid_from_name(extract_species_name(x) or "", taxdb)[0])
    ictv["IDS"] = ictv.apply(
        lambda x: extract_from_hyperlink(x["Accessions Link"])["url"].split("/")[-1].split(","),
        axis=1,
    )
    ictv["Range"] = ictv.apply(lambda x: extract_regions(x["Virus GENBANK accession"]), axis=1)

    id2newtax = {}
    id2range = {}
    for _, row in ictv[["Taxid", "Virus name(s)", "IDS", "Range"]].iterrows():
        for acc in row["IDS"]:
            acc = acc.split("(")[0]
            id2newtax[acc] = row["Taxid"]
            if row["Range"] and len(row["Range"]) > 0:
                id2range[acc] = row["Range"]

    ictv["IDS_NR"] = ictv["IDS"].apply(filter_genbank)

    all_nrids = []
    for _, row in ictv.iterrows():
        if len(row["IDS_NR"]) > 0 and row["IDS_NR"][0] != "":
            all_nrids.extend(row["IDS_NR"])

    all_nrids = set(all_nrids)
    all_files = glob(str(gbdir / "*.gb"))

    proteins = []
    for file in tqdm(all_files, desc="Extracting proteins"):
        with open(file) as handle:
            for record in GenBank.parse(handle):
                if record.version.split(".")[0] not in all_nrids:  # type: ignore
                    continue

                if record.version.split(".")[0] in id2range:  # type: ignore
                    feature_range = id2range[record.version.split(".")[0]][0]  # type: ignore
                else:
                    feature_range = (0, np.inf)

                for feature in record.features:  # type: ignore
                    if feature.key != "CDS":
                        continue

                    protein_id = ""
                    product = ""
                    protein_seq = ""

                    for qualifier in feature.qualifiers:
                        if qualifier.key == "/protein_id=":
                            protein_id = qualifier.value.strip('"').rstrip('"')
                        elif qualifier.key == "/product=":
                            product = "[" + qualifier.value.strip('"').rstrip('"') + "]"
                        elif qualifier.key == "/translation=":
                            protein_seq = qualifier.value.strip('"').rstrip('"')

                    start, end = list(map(int, re.findall(r"(\d+)\D+(\d+)", feature.location)[0]))
                    if (min(feature_range[0], start) == feature_range[0] and max(feature_range[1], end) == feature_range[1]) and len(protein_seq) > 0:
                        proteins.append(
                            SeqRecord(
                                Seq(protein_seq),
                                id=protein_id,
                                description=product,
                                name=product,
                            )
                        )

    proteins_path = seqdir / "proteins.faa"
    with open(proteins_path, "w") as out:
        SeqIO.write(proteins, out, "fasta")
    logger.info("%s created", proteins_path)

    prot_map_path = seqdir / "virus_protein.accession2taxid"
    with open(prot_map_path, "w") as fh2:
        fh2.write("accession\taccession.version\ttaxid\tgi\n")
        for file in tqdm(all_files, desc="Writing protein accession2taxid"):
            with open(file) as handle:
                for record in GenBank.parse(handle):
                    if record.version.split(".")[0] not in all_nrids:  # type: ignore
                        continue

                    if record.version.split(".")[0] in id2range:  # type: ignore
                        feature_range = id2range[record.version.split(".")[0]][0]  # type: ignore
                    else:
                        feature_range = (0, np.inf)

                    for feature in record.features:  # type: ignore
                        if feature.key != "CDS":
                            continue

                        protein_id = ""
                        protein_seq = ""

                        for qualifier in feature.qualifiers:
                            if qualifier.key == "/protein_id=":
                                protein_id = qualifier.value.strip('"').rstrip('"')
                            elif qualifier.key == "/translation=":
                                protein_seq = qualifier.value.strip('"').rstrip('"')

                        start, end = list(map(int, re.findall(r"(\d+)\D+(\d+)", feature.location)[0]))
                        if (min(feature_range[0], start) == feature_range[0] and max(feature_range[1], end) == feature_range[1]) and len(protein_seq) > 0:
                            fh2.write(
                                f"{protein_id.split('.')[0]}\t{protein_id}\t{id2newtax[record.accession[0]]}\t-\n"  # type: ignore
                            )
    logger.info("%s created", prot_map_path)

    seq_records = {}
    for file in tqdm(all_files, desc="Extracting genomes"):
        with open(file) as handle:
            for record in GenBank.parse(handle):
                if record.version.split(".")[0] not in all_nrids:  # type: ignore
                    continue

                if record.version.split(".")[0] in id2range:  # type: ignore
                    feature_range = id2range[record.version.split(".")[0]][0]  # type: ignore
                else:
                    feature_range = (0, None)

                seq = record.sequence[feature_range[0] : feature_range[1]]  # type: ignore
                if not seq:
                    continue

                if record.accession[0] in seq_records:  # type: ignore
                    logger.warning("%s is duplicated!", record.accession[0])  # type: ignore
                    continue

                seq_records[record.accession[0]] = SeqRecord(  # type: ignore
                    Seq(seq),
                    id=record.accession[0],  # type: ignore
                    description=f"[{record.source}] {record.data_file_division}",  # type: ignore
                )

    genomes_path = seqdir / "genomes.fna"
    with open(genomes_path, "w") as out:
        SeqIO.write(seq_records.values(), out, "fasta")
    logger.info("%s created", genomes_path)

    genome_map_path = seqdir / "virus_genome.accession2taxid"
    with open(genome_map_path, "w") as fh2:
        fh2.write("accession\taccession.version\ttaxid\tgi\n")
        for file in tqdm(all_files, desc="Writing genome accession2taxid"):
            with open(file) as handle:
                for record in GenBank.parse(handle):
                    if record.version.split(".")[0] in all_nrids: # type: ignore
                        fh2.write(
                            f"{record.accession[0].split('.')[0]}\t{record.accession[0].split('.')[0]}\t{id2newtax[record.accession[0]]}\t-\n"  # type: ignore
                        )

    logger.info("%s created", genome_map_path)


if __name__ == "__main__":
    main()
