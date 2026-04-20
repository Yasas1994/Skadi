#!/usr/bin/env python3

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import click
from Bio import SeqIO


DEFAULT_EXTENSIONS = {".gb", ".gbk", ".gbff", ".genbank"}


def discover_genbank_files(inputs: list[Path], exts: set[str]) -> list[Path]:
    files = []

    for p in inputs:
        if p.is_file():
            if p.suffix.lower() in exts:
                files.append(p)
        elif p.is_dir():
            for child in sorted(p.rglob("*")):
                if child.is_file() and child.suffix.lower() in exts:
                    files.append(child)

    # unique + stable order
    seen = set()
    out = []
    for f in files:
        rf = f.resolve()
        if rf not in seen:
            seen.add(rf)
            out.append(f)

    return out


def first_qualifier(feature, key: str, default: str = "") -> str:
    value = feature.qualifiers.get(key, [])
    if not value:
        return default
    if isinstance(value, list):
        return str(value[0])
    return str(value)


def get_record_id(record) -> str:
    accessions = record.annotations.get("accessions", [])
    if accessions:
        return str(accessions[0])
    if getattr(record, "id", None):
        return str(record.id)
    if getattr(record, "name", None):
        return str(record.name)
    return "unknown_record"


def extract_rows_from_file(gb_file: Path, group_by: str) -> list[dict]:
    rows = []

    for record in SeqIO.parse(str(gb_file), "genbank"):
        record_id = get_record_id(record)
        #genome_id = gb_file.stem if group_by == "file" else record_id

        for feature in record.features:
            if feature.type != "CDS":
                continue

            protein_ids = feature.qualifiers.get("protein_id", [])
            if not protein_ids:
                continue

            locus_tag = first_qualifier(feature, "locus_tag")
            gene = first_qualifier(feature, "gene")
            product = first_qualifier(feature, "product")

            for protein_id in protein_ids:
                rows.append(
                    {
                        "genome_id": record_id,
                        "protein_id": str(protein_id),
                        "locus_tag": locus_tag,
                        "gene": gene,
                        "product": product,
                    }
                )

    return rows


def deduplicate_rows(rows: list[dict]) -> list[dict]:
    seen = set()
    out = []

    for row in rows:
        key = (
            row["genome_id"],
            row["protein_id"],
            row["locus_tag"],
        )
        if key not in seen:
            seen.add(key)
            out.append(row)

    return out


def write_tsv(rows: list[dict], output_file: Path) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)

    header = ["genome_id", "protein_id", "locus_tag", "gene", "product"]

    with output_file.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(row.get(col, "") for col in header) + "\n")


def build_json_map(rows: list[dict]) -> dict[str, list[str]]:
    mapping = defaultdict(set)
    for row in rows:
        mapping[row["genome_id"]].add(row["protein_id"])

    return {genome_id: sorted(protein_ids) for genome_id, protein_ids in sorted(mapping.items())}


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument(
    "inputs",
    nargs=-1,
    required=True,
    type=click.Path(exists=True, path_type=Path),
)
@click.option(
    "-o",
    "--output",
    "output_file",
    required=True,
    type=click.Path(path_type=Path),
    help="Output TSV file.",
)
@click.option(
    "--json-output",
    type=click.Path(path_type=Path),
    default=None,
    help="Optional JSON output of genome_id -> [protein_id, ...].",
)
@click.option(
    "--group-by",
    type=click.Choice(["file", "record"]),
    default="file",
    show_default=True,
    help="How to define genome_id.",
)
@click.option(
    "--extensions",
    default=".gb,.gbk,.gbff,.genbank",
    show_default=True,
    help="Comma-separated file extensions to scan.",
)
@click.option(
    "--deduplicate/--keep-duplicates",
    default=True,
    show_default=True,
    help="Deduplicate repeated genome_id/protein_id entries.",
)
def cli(
    inputs: tuple[Path, ...],
    output_file: Path,
    json_output: Path | None,
    group_by: str,
    extensions: str,
    deduplicate: bool,
) -> None:
    """
    Create a genome_id -> protein_id map from one or more GenBank files.

    INPUTS can be one or more GenBank files and/or directories.
    """
    exts = {x.strip().lower() for x in extensions.split(",") if x.strip()}
    gb_files = discover_genbank_files(list(inputs), exts)

    if not gb_files:
        raise click.ClickException("No GenBank files found.")

    all_rows = []
    for gb_file in gb_files:
        rows = extract_rows_from_file(gb_file, group_by=group_by)
        all_rows.extend(rows)

    if deduplicate:
        all_rows = deduplicate_rows(all_rows)

    if not all_rows:
        raise click.ClickException("No CDS features with protein_id qualifiers were found.")

    write_tsv(all_rows, output_file)
    click.echo(f"Wrote TSV: {output_file} ({len(all_rows)} rows)")

    if json_output is not None:
        mapping = build_json_map(all_rows)
        json_output.parent.mkdir(parents=True, exist_ok=True)
        with json_output.open("w", encoding="utf-8") as fh:
            json.dump(mapping, fh, indent=2)
        click.echo(f"Wrote JSON: {json_output} ({len(mapping)} genome IDs)")


if __name__ == "__main__":
    cli()