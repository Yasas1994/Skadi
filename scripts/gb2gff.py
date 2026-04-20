#!/usr/bin/env python3
"""
Convert a GenBank file that may contain multiple records into
one GFF3 file per record, using the primary accession as the filename.

It verifies before exiting that every GenBank record accession has a
matching GFF3 file in the output directory.

Dependencies:
    pip install click biopython
"""

from __future__ import annotations

import re
from pathlib import Path
from urllib.parse import quote

import click
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation


def sanitize_filename(name: str) -> str:
    """Make a safe filename stem."""
    name = (name or "").strip() or "record"
    name = re.sub(r"[^\w.\-]+", "_", name)
    return name[:200]


def gff3_escape(value) -> str:
    """Escape a value for GFF3 attributes."""
    return quote(str(value), safe="._:-|")


def flatten_qualifier_value(value) -> str:
    """Convert GenBank qualifier values into a GFF3-safe string."""
    if isinstance(value, (list, tuple)):
        return ",".join(gff3_escape(v) for v in value)
    return gff3_escape(value)


def get_primary_accession(record) -> str:
    """
    Return the primary accession for a GenBank record.

    Priority:
      1. record.annotations['accessions'][0]
      2. record.id
      3. record.name
    """
    accessions = record.annotations.get("accessions", [])
    if accessions:
        acc = str(accessions[0]).strip()
        if acc and acc not in {"<unknown id>", "<unknown name>"}:
            return sanitize_filename(acc)

    for candidate in (getattr(record, "id", None), getattr(record, "name", None)):
        if candidate:
            candidate = str(candidate).strip()
            if candidate and candidate not in {"<unknown id>", "<unknown name>"}:
                # record.id can sometimes include version, which is usually fine
                return sanitize_filename(candidate)

    return "record"


def output_path_for_accession(output_dir: Path, accession: str) -> Path:
    """Return the exact expected GFF3 path for an accession."""
    return output_dir / f"{sanitize_filename(accession)}.gff3"


def get_source(feature) -> str:
    source_values = feature.qualifiers.get("source")
    if source_values:
        if isinstance(source_values, list) and source_values:
            return str(source_values[0])
        return str(source_values)
    return "GenBank"


def feature_name(feature) -> str | None:
    for key in ("gene", "locus_tag", "product", "protein_id", "note"):
        value = feature.qualifiers.get(key)
        if value:
            if isinstance(value, list):
                return str(value[0])
            return str(value)
    return None


def build_attributes(feature, feature_id: str, record_id: str, accession: str) -> str:
    attrs = {"ID": feature_id}

    name = feature_name(feature)
    if name:
        attrs["Name"] = name

    for key, value in feature.qualifiers.items():
        if key in {"translation"}:
            continue
        if key == "db_xref":
            attrs["Dbxref"] = flatten_qualifier_value(value)
        elif key == "note":
            attrs["Note"] = flatten_qualifier_value(value)
        elif key == "gene":
            attrs["gene"] = flatten_qualifier_value(value)
        elif key == "locus_tag":
            attrs["locus_tag"] = flatten_qualifier_value(value)
        elif key == "product":
            attrs["product"] = flatten_qualifier_value(value)
        elif key == "protein_id":
            attrs["protein_id"] = flatten_qualifier_value(value)
        elif key == "gene_synonym":
            attrs["Alias"] = flatten_qualifier_value(value)
        elif key == "codon_start":
            attrs["codon_start"] = flatten_qualifier_value(value)

    if feature.type == "source":
        if "Name" not in attrs:
            attrs["Name"] = gff3_escape(record_id)
        attrs["accession"] = gff3_escape(accession)

    return ";".join(f"{k}={v}" for k, v in attrs.items())


def location_parts(feature):
    loc = feature.location
    if isinstance(loc, CompoundLocation):
        return list(loc.parts)
    return [loc]


def get_phase_map_for_cds(feature) -> dict[tuple[int, int, int | None], str]:
    """
    Compute GFF3 phase for each CDS segment.

    GFF3 phase:
      0 -> segment starts at codon boundary
      1 -> skip first base to reach codon boundary
      2 -> skip first 2 bases to reach codon boundary
    """
    parts = location_parts(feature)
    if not parts:
        return {}

    genomic_parts = sorted(parts, key=lambda p: int(p.start))
    strand = getattr(feature.location, "strand", None)

    if strand == -1:
        tx_parts = list(reversed(genomic_parts))
    else:
        tx_parts = genomic_parts

    codon_start = 1
    raw_codon_start = feature.qualifiers.get("codon_start")
    if raw_codon_start:
        try:
            codon_start = int(raw_codon_start[0] if isinstance(raw_codon_start, list) else raw_codon_start)
            if codon_start not in (1, 2, 3):
                codon_start = 1
        except Exception:
            codon_start = 1

    phase_map = {}
    consumed = 0

    for i, part in enumerate(tx_parts):
        part_len = int(part.end) - int(part.start)

        if i == 0:
            phase = codon_start - 1
            coding_len = max(part_len - phase, 0)
        else:
            phase = (3 - (consumed % 3)) % 3
            coding_len = max(part_len - phase, 0)

        key = (int(part.start), int(part.end), part.strand)
        phase_map[key] = str(phase)
        consumed += coding_len

    return phase_map


def strand_symbol(strand) -> str:
    if strand == 1:
        return "+"
    if strand == -1:
        return "-"
    return "."


def write_record_gff3(record, output_path: Path, accession: str) -> None:
    seqid = record.id
    seq_len = len(record.seq)

    with output_path.open("w", encoding="utf-8") as out:
        out.write("##gff-version 3\n")
        out.write(f"##sequence-region {gff3_escape(seqid)} 1 {seq_len}\n")

        feature_counts: dict[str, int] = {}

        for feature in record.features:
            ftype = feature.type or "feature"
            feature_counts[ftype] = feature_counts.get(ftype, 0) + 1
            feature_id = f"{ftype}_{feature_counts[ftype]}"

            source = gff3_escape(get_source(feature))
            ftype_escaped = gff3_escape(ftype)
            score = "."
            attrs = build_attributes(feature, feature_id, seqid, accession)

            parts = sorted(location_parts(feature), key=lambda p: int(p.start))
            phase_map = get_phase_map_for_cds(feature) if ftype == "CDS" else {}

            for part in parts:
                start = int(part.start) + 1
                end = int(part.end)
                key = (int(part.start), int(part.end), part.strand)
                phase = phase_map.get(key, ".") if ftype == "CDS" else "."

                out.write(
                    "\t".join(
                        [
                            gff3_escape(seqid),
                            source,
                            ftype_escaped,
                            str(start),
                            str(end),
                            score,
                            strand_symbol(part.strand),
                            phase,
                            attrs,
                        ]
                    )
                    + "\n"
                )


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument(
    "genbank_file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.argument(
    "output_dir",
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite existing GFF3 files if they already exist.",
)
def cli(genbank_file: Path, output_dir: Path, overwrite: bool) -> None:
    """
    Convert each GenBank record in GENBANK_FILE into a separate GFF3 file
    inside OUTPUT_DIR, naming files by accession and verifying all expected
    accession-matched GFF3 files exist before exiting.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    records = list(SeqIO.parse(str(genbank_file), "genbank"))
    if not records:
        raise click.ClickException("No GenBank records were found in the input file.")

    expected_accessions: list[str] = []
    seen_accessions: set[str] = set()

    for record in records:
        accession = get_primary_accession(record)

        if accession in seen_accessions:
            click.echo(
                f"Duplicate accession detected in input GenBank file: {accession}"
            )

        seen_accessions.add(accession)
        expected_accessions.append(accession)

        out_path = output_path_for_accession(output_dir, accession)

        if out_path.exists() and not overwrite:
            click.echo(f"Exists, skipping: {out_path}")
            continue

        write_record_gff3(record, out_path, accession)
        click.echo(f"Wrote: {out_path}")

    missing = []
    for accession in expected_accessions:
        expected_path = output_path_for_accession(output_dir, accession)
        if not expected_path.exists():
            missing.append((accession, expected_path))

    if missing:
        msg = "\n".join(f"{acc}\t{path}" for acc, path in missing[:20])
        extra = ""
        if len(missing) > 20:
            extra = f"\n... and {len(missing) - 20} more"
        raise click.ClickException(
            "Verification failed: some expected accession-matched GFF3 files are missing:\n"
            f"{msg}{extra}"
        )

    click.echo(
        f"Done. Verified {len(expected_accessions)} accession-matched GFF3 file(s) in {output_dir}"
    )


if __name__ == "__main__":
    cli()