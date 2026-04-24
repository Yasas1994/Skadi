#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import click


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def guess_delimiter(line: str) -> str:
    if "\t" in line:
        return "\t"
    if "," in line:
        return ","
    return "\t"


def normalize_token(x: str) -> str:
    return x.strip().split()[0]


def strip_version(x: str) -> str:
    return re.sub(r"\.\d+$", "", x)


def build_alias_lookup(keys: Iterable[str]) -> dict[str, str]:
    alias_to_real: dict[str, str] = {}
    ambiguous: set[str] = set()

    for key in keys:
        token = normalize_token(key)
        aliases = {token, strip_version(token)}

        for alias in aliases:
            if not alias or alias == key:
                continue
            if alias in ambiguous:
                continue
            if alias in alias_to_real and alias_to_real[alias] != key:
                ambiguous.add(alias)
                alias_to_real.pop(alias, None)
            else:
                alias_to_real[alias] = key

    return alias_to_real


def resolve_id(x: str, exact: dict, alias: dict):
    if x in exact:
        return x

    token = normalize_token(x)
    if token in exact:
        return token

    if x in alias:
        return alias[x]
    if token in alias:
        return alias[token]

    sv = strip_version(token)
    if sv in exact:
        return sv
    if sv in alias:
        return alias[sv]

    return None


def parse_attributes(attr_str: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    if not attr_str or attr_str == ".":
        return attrs

    for item in attr_str.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k] = v
        else:
            attrs[item] = "true"
    return attrs


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []

    intervals = sorted((min(a, b), max(a, b)) for a, b in intervals)
    merged = [intervals[0]]

    for start, end in intervals[1:]:
        cur_start, cur_end = merged[-1]
        if start <= cur_end + 1:
            merged[-1] = (cur_start, max(cur_end, end))
        else:
            merged.append((start, end))

    return merged


def union_length(intervals: list[tuple[int, int]]) -> int:
    return sum(end - start + 1 for start, end in merge_intervals(intervals))


def sum_length(intervals: list[tuple[int, int]]) -> int:
    return sum(max(0, end - start + 1) for start, end in intervals)


def safe_div(num: float, den: float) -> float:
    return num / den if den else 0.0


RANK_ORDER = [
    "r", "sr", "k", "sk", "p", "sp", "c", "sc", "o", "so", "f", "sf", "g", "sg", "s"
]

PREFIX_RE = re.compile(r"^(r|sr|k|sk|p|sp|c|sc|o|so|f|sf|g|sg|s)__", re.IGNORECASE)

SUFFIX_TO_RANK = [
    ("viricotina", "sp"),
    ("viricetidae", "sc"),
    ("viricetes", "c"),
    ("viricota", "p"),
    ("virineae", "so"),
    ("virinae", "sf"),
    ("virites", "sk"),
    ("virales", "o"),
    ("viridae", "f"),
    ("virae", "k"),
    ("viria", "r"),
    ("vira", "sr"),
    ("virus", "g"),
]


def clean_taxon_name(token: str) -> str:
    token = token.strip().strip(";")
    if not token:
        return ""

    if PREFIX_RE.match(token):
        token = token.split("__", 1)[1]
    elif "_" in token:
        token = token.split("_", 1)[1]

    token = token.strip().strip(";")
    token = token.replace("_", " ")
    token = re.sub(r"\s+", " ", token)
    return token.strip()


def classify_viral_token(name: str, ranks: dict[str, str]) -> None:
    lname = name.lower()

    for suffix, rank in SUFFIX_TO_RANK:
        if lname.endswith(suffix):
            if rank == "g":
                if "g" not in ranks:
                    ranks["g"] = name
                elif "sg" not in ranks and ranks["g"] != name:
                    ranks["sg"] = name
            else:
                if rank not in ranks:
                    ranks[rank] = name
            return

    if " " in name and "s" not in ranks:
        ranks["s"] = name


def format_taxlineage(raw_lineage: str) -> str:
    if raw_lineage is None:
        return ""

    raw_lineage = raw_lineage.strip()
    if raw_lineage == "":
        return ""

    tokens = [t.strip() for t in raw_lineage.split(";") if t.strip()]
    if not tokens:
        return ""

    ranks: dict[str, str] = {}

    has_prefixed = any(PREFIX_RE.match(t) for t in tokens)
    if has_prefixed:
        for token in tokens:
            m = PREFIX_RE.match(token)
            if not m:
                continue
            prefix = m.group(1).lower()
            value = clean_taxon_name(token)
            if value and prefix not in ranks:
                ranks[prefix] = value
    else:
        for token in tokens:
            name = clean_taxon_name(token)
            if not name:
                continue
            classify_viral_token(name, ranks)

    return ";".join(f"{rank}__{ranks[rank]}" for rank in RANK_ORDER if rank in ranks and ranks[rank])


def normalize_staxid(raw: str) -> str:
    raw = (raw or "").strip()
    if raw == "":
        return ""
    for token in re.split(r"[;,]", raw):
        token = token.strip()
        if token:
            return token
    return ""


def read_protein_map(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    first_line = None
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                first_line = line
                break

    if first_line is None:
        raise click.ClickException(f"Protein map is empty: {path}")

    delim = guess_delimiter(first_line)
    header_fields = [x.strip().lower() for x in first_line.split(delim)]

    protein_to_genome: dict[str, str] = {}
    has_header = ("protein_id" in header_fields and "genome_id" in header_fields)

    with open_text(path) as fh:
        reader = csv.reader(fh, delimiter=delim)

        if has_header:
            header = None
            for row in reader:
                if not row or row[0].startswith("#"):
                    continue
                header = [x.strip().lower() for x in row]
                break

            assert header is not None
            pidx = header.index("protein_id")
            gidx = header.index("genome_id")

            for row in reader:
                if not row or row[0].startswith("#"):
                    continue
                if len(row) <= max(pidx, gidx):
                    continue

                protein_id = row[pidx].strip()
                genome_id = row[gidx].strip()

                if protein_id and genome_id:
                    protein_to_genome[protein_id] = genome_id
        else:
            for row in reader:
                if not row or row[0].startswith("#"):
                    continue
                if len(row) < 2:
                    continue

                protein_id = row[0].strip()
                genome_id = row[1].strip()

                if protein_id and genome_id:
                    protein_to_genome[protein_id] = genome_id

    alias = build_alias_lookup(protein_to_genome.keys())
    return protein_to_genome, alias


@dataclass
class ProteinModel:
    protein_id: str
    seqid: str
    strand: str
    segments: list[tuple[int, int]] = field(default_factory=list)

    @property
    def tx_segments(self) -> list[tuple[int, int]]:
        if self.strand == "-":
            return sorted(self.segments, key=lambda x: x[0], reverse=True)
        return sorted(self.segments, key=lambda x: x[0])

    @property
    def cds_nt_len(self) -> int:
        return sum(end - start + 1 for start, end in self.segments)

    @property
    def aa_len(self) -> int:
        return self.cds_nt_len // 3

    def aa_interval_to_genome_intervals(self, aa_start: int, aa_end: int) -> list[tuple[int, int]]:
        if self.aa_len <= 0:
            return []

        aa_start = max(1, aa_start)
        aa_end = min(self.aa_len, aa_end)
        if aa_start > aa_end:
            return []

        nt_start0 = (aa_start - 1) * 3
        nt_end0 = aa_end * 3

        out: list[tuple[int, int]] = []
        tx_pos0 = 0

        for seg_start, seg_end in self.tx_segments:
            seg_len = seg_end - seg_start + 1
            seg_tx_start0 = tx_pos0
            seg_tx_end0 = tx_pos0 + seg_len

            ov_start0 = max(nt_start0, seg_tx_start0)
            ov_end0 = min(nt_end0, seg_tx_end0)

            if ov_start0 < ov_end0:
                local_start0 = ov_start0 - seg_tx_start0
                local_end0 = ov_end0 - seg_tx_start0

                if self.strand == "-":
                    g_start = seg_end - local_end0 + 1
                    g_end = seg_end - local_start0
                else:
                    g_start = seg_start + local_start0
                    g_end = seg_start + local_end0 - 1

                out.append((min(g_start, g_end), max(g_start, g_end)))

            tx_pos0 += seg_len

        return out


def parse_gff3_models(
    gff3_path: Path,
    protein_attr_keys: list[str],
) -> tuple[dict[str, ProteinModel], dict[str, str], list[tuple[int, int]]]:
    proteins: dict[str, ProteinModel] = {}
    all_cds_intervals: list[tuple[int, int]] = []

    with open_text(gff3_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            if ftype != "CDS":
                continue

            start_i = int(start)
            end_i = int(end)
            all_cds_intervals.append((min(start_i, end_i), max(start_i, end_i)))

            attr_d = parse_attributes(attrs)

            protein_id = None
            for key in protein_attr_keys:
                if key in attr_d and attr_d[key]:
                    protein_id = attr_d[key].split(",")[0]
                    break

            if protein_id is None:
                continue

            if protein_id not in proteins:
                proteins[protein_id] = ProteinModel(
                    protein_id=protein_id,
                    seqid=seqid,
                    strand=strand,
                    segments=[(min(start_i, end_i), max(start_i, end_i))],
                )
            else:
                proteins[protein_id].segments.append((min(start_i, end_i), max(start_i, end_i)))

    alias = build_alias_lookup(proteins.keys())
    return proteins, alias, all_cds_intervals


def read_blast_hits(
    blast_path: Path,
    blast_columns: list[str],
    min_pident: float | None,
    max_evalue: float | None,
    min_bitscore: float | None,
):
    required = {
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "staxids",
        "slineages",
    }
    if not required.issubset(set(blast_columns)):
        missing = sorted(required - set(blast_columns))
        raise click.ClickException(
            f"--blast-columns is missing required fields: {', '.join(missing)}"
        )

    first_line = None
    with open_text(blast_path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                first_line = line
                break

    if first_line is None:
        return

    delim = "\t"
    first_fields = [x.strip() for x in first_line.split(delim)]
    has_header = (
        set(first_fields) >= required
        or "qseqid" in first_fields
        or "sseqid" in first_fields
    )

    with open_text(blast_path) as fh:
        reader = csv.reader(fh, delimiter=delim)
        header = blast_columns

        if has_header:
            for row in reader:
                if not row or row[0].startswith("#"):
                    continue
                header = [x.strip() for x in row]
                break

        h2i = {name: i for i, name in enumerate(header)}

        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) < len(header):
                continue

            rec: dict[str, str | int | float] = {name: row[idx] for name, idx in h2i.items() if idx < len(row)}

            try:
                rec["pident"] = float(rec["pident"])
                rec["length"] = int(rec["length"])
                rec["sstart"] = int(rec["sstart"])
                rec["send"] = int(rec["send"])
                rec["evalue"] = float(rec["evalue"])
                rec["bitscore"] = float(rec["bitscore"])
            except ValueError:
                continue

            if min_pident is not None and rec["pident"] < min_pident:
                continue

            if max_evalue is not None and rec["evalue"] > max_evalue:
                continue

            if min_bitscore is not None and rec["bitscore"] < min_bitscore:
                continue

            yield rec


def index_gff3_files(gff3_dir: Path, extension: str) -> dict[str, Path]:
    files = sorted(gff3_dir.rglob(f"*{extension}"))
    index: dict[str, Path] = {}

    for path in files:
        stem = path.stem
        if stem in index:
            raise click.ClickException(
                f"Duplicate GFF3 stem '{stem}' found:\n  {index[stem]}\n  {path}\n"
                "Genome IDs are matched to GFF3 file stems."
            )
        index[stem] = path

    return index


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("blast_tsv", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.argument("protein_map", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.argument("gff3_dir", type=click.Path(exists=True, file_okay=False, path_type=Path))
@click.argument("output_tsv", type=click.Path(path_type=Path))
@click.option(
    "--blast-columns",
    default="qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,staxids,slineages",
    show_default=True,
    help="Comma-separated DIAMOND/BLAST outfmt 6 columns used if the blast file has no header.",
)
@click.option(
    "--protein-attr",
    default="protein_id,Name",
    show_default=True,
    help="Comma-separated GFF3 attribute keys to try for protein IDs.",
)
@click.option(
    "--gff3-ext",
    default=".gff3",
    show_default=True,
    help="GFF3 file extension to scan for under gff3_dir.",
)
@click.option("--min-pident", type=float, default=None, help="Minimum percent identity filter.")
@click.option("--max-evalue", type=float, default=None, help="Maximum e-value filter.")
@click.option("--min-bitscore", type=float, default=None, help="Minimum bitscore filter.")
@click.option(
    "--min-virus-protein-coverage",
    type=float,
    default=None,
    help="Minimum proteins_hit / proteins_total required to keep a genome in the output.",
)
@click.option(
    "--min-avg-aa-identity-weighted",
    type=float,
    default=None,
    help="Minimum avg_aa_identity_weighted required to keep a genome in the output.",
)
@click.option(
    "--include-zero-hit-genomes/--only-hit-genomes",
    default=False,
    show_default=True,
    help="Include genomes with GFF3 files but zero retained hits.",
)
def cli(
    blast_tsv: Path,
    protein_map: Path,
    gff3_dir: Path,
    output_tsv: Path,
    blast_columns: str,
    protein_attr: str,
    gff3_ext: str,
    min_pident: float | None,
    max_evalue: float | None,
    min_bitscore: float | None,
    min_virus_protein_coverage: float | None,
    min_avg_aa_identity_weighted: float | None,
    include_zero_hit_genomes: bool,
):
    blast_cols = [x.strip() for x in blast_columns.split(",") if x.strip()]
    protein_attr_keys = [x.strip() for x in protein_attr.split(",") if x.strip()]

    if min_virus_protein_coverage is not None:
        if not (0.0 <= min_virus_protein_coverage <= 1.0):
            raise click.ClickException("--min-virus-protein-coverage must be between 0 and 1.")

    protein_to_genome, protein_map_alias = read_protein_map(protein_map)
    gff_index = index_gff3_files(gff3_dir, gff3_ext)

    hits_by_genome: dict[str, list[dict]] = defaultdict(list)
    taxonomy_by_genome: dict[str, Counter] = defaultdict(Counter)

    blast_rows_total = 0
    blast_rows_kept = 0
    blast_rows_unmapped_subject = 0

    protein_exact_lookup = {k: True for k in protein_to_genome.keys()}

    for rec in read_blast_hits(
        blast_path=blast_tsv,
        blast_columns=blast_cols,
        min_pident=min_pident,
        max_evalue=max_evalue,
        min_bitscore=min_bitscore,
    ) or []:
        blast_rows_total += 1

        resolved_protein = resolve_id(str(rec["sseqid"]), protein_exact_lookup, protein_map_alias)
        if resolved_protein is None:
            blast_rows_unmapped_subject += 1
            continue

        genome_id = protein_to_genome[resolved_protein]
        aa_start = min(rec["sstart"], rec["send"])
        aa_end = max(rec["sstart"], rec["send"])

        taxid = normalize_staxid(str(rec.get("staxids", "")))
        taxlineage = format_taxlineage(str(rec.get("slineages", "")))

        hits_by_genome[genome_id].append(
            {
                "qseqid": rec["qseqid"],
                "protein_id": resolved_protein,
                "aa_start": aa_start,
                "aa_end": aa_end,
                "taxid": taxid,
                "taxlineage": taxlineage,
                "pident": rec["pident"],
                "align_len_aa": rec["length"],
            }
        )

        taxonomy_by_genome[genome_id][(taxid, taxlineage)] += 1
        blast_rows_kept += 1

    genomes_to_process = set(hits_by_genome.keys())
    if include_zero_hit_genomes:
        genomes_to_process.update(gff_index.keys())

    if not genomes_to_process:
        raise click.ClickException("No genomes to process after filtering and mapping.")

    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "genome_id",
        "taxid",
        "taxlineage",
        "proteins_total",
        "proteins_hit",
        "virus_protein_coverage",
        "queries_hit",
        "blast_hits",
        "coding_nt_total",
        "coding_nt_covered",
        "genome_coverage",
        "genome_mean_depth",
        "genome_mean_depth_on_covered",
        "proteome_aa_total",
        "proteome_aa_covered",
        "protein_coverage",
        "protein_mean_depth",
        "protein_mean_depth_on_covered",
        "avg_aa_identity",
        "avg_aa_identity_weighted",
        "hits_missing_in_gff",
    ]

    rows_out: list[dict] = []
    missing_gff = 0
    missing_protein_in_gff = 0
    filtered_by_protein_coverage = 0
    filtered_by_weighted_identity = 0

    for genome_id in sorted(genomes_to_process):
        gff3_path = gff_index.get(genome_id)
        if gff3_path is None:
            missing_gff += 1
            click.echo(f"Warning: no GFF3 found for genome_id={genome_id}", err=True)
            continue

        taxid = ""
        taxlineage = ""
        if taxonomy_by_genome.get(genome_id):
            (taxid, taxlineage), _ = taxonomy_by_genome[genome_id].most_common(1)[0]

        proteins, proteins_alias, all_cds_intervals = parse_gff3_models(
            gff3_path,
            protein_attr_keys=protein_attr_keys,
        )

        protein_exact_lookup_this = {k: True for k in proteins.keys()}

        proteome_aa_total = sum(p.aa_len for p in proteins.values())
        coding_nt_total = union_length(all_cds_intervals)

        protein_hit_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
        genome_hit_intervals: list[tuple[int, int]] = []
        queries_hit: set[str] = set()
        hits_missing_here = 0

        genome_hits = hits_by_genome.get(genome_id, [])

        avg_aa_identity = 0.0
        avg_aa_identity_weighted = 0.0
        if genome_hits:
            avg_aa_identity = sum(hit["pident"] for hit in genome_hits) / len(genome_hits)
            weighted_den = sum(max(0, hit["align_len_aa"]) for hit in genome_hits)
            if weighted_den > 0:
                avg_aa_identity_weighted = (
                    sum(hit["pident"] * max(0, hit["align_len_aa"]) for hit in genome_hits)
                    / weighted_den
                )

        for hit in genome_hits:
            queries_hit.add(hit["qseqid"])

            resolved_protein = resolve_id(
                hit["protein_id"],
                protein_exact_lookup_this,
                proteins_alias,
            )
            if resolved_protein is None:
                hits_missing_here += 1
                missing_protein_in_gff += 1
                continue

            model = proteins[resolved_protein]
            if model.aa_len <= 0:
                continue

            aa_start = max(1, hit["aa_start"])
            aa_end = min(model.aa_len, hit["aa_end"])
            if aa_start > aa_end:
                continue

            protein_hit_intervals[resolved_protein].append((aa_start, aa_end))
            genome_hit_intervals.extend(model.aa_interval_to_genome_intervals(aa_start, aa_end))

        proteome_aa_covered = sum(union_length(v) for v in protein_hit_intervals.values())
        proteome_aa_depth_sum = sum(sum_length(v) for v in protein_hit_intervals.values())

        coding_nt_covered = union_length(genome_hit_intervals)
        coding_nt_depth_sum = sum_length(genome_hit_intervals)

        proteins_total = len(proteins)
        proteins_hit = len(protein_hit_intervals)
        virus_protein_coverage = safe_div(proteins_hit, proteins_total)

        if min_virus_protein_coverage is not None and virus_protein_coverage < min_virus_protein_coverage:
            filtered_by_protein_coverage += 1
            continue

        if (
            min_avg_aa_identity_weighted is not None
            and avg_aa_identity_weighted < min_avg_aa_identity_weighted
        ):
            filtered_by_weighted_identity += 1
            continue

        rows_out.append(
            {
                "genome_id": genome_id,
                "taxid": taxid,
                "taxlineage": taxlineage,
                "proteins_total": proteins_total,
                "proteins_hit": proteins_hit,
                "virus_protein_coverage": f"{virus_protein_coverage:.6f}",
                "queries_hit": len(queries_hit),
                "blast_hits": len(genome_hits),
                "coding_nt_total": coding_nt_total,
                "coding_nt_covered": coding_nt_covered,
                "genome_coverage": f"{safe_div(coding_nt_covered, coding_nt_total):.6f}",
                "genome_mean_depth": f"{safe_div(coding_nt_depth_sum, coding_nt_total):.6f}",
                "genome_mean_depth_on_covered": f"{safe_div(coding_nt_depth_sum, coding_nt_covered):.6f}",
                "proteome_aa_total": proteome_aa_total,
                "proteome_aa_covered": proteome_aa_covered,
                "protein_coverage": f"{safe_div(proteome_aa_covered, proteome_aa_total):.6f}",
                "protein_mean_depth": f"{safe_div(proteome_aa_depth_sum, proteome_aa_total):.6f}",
                "protein_mean_depth_on_covered": f"{safe_div(proteome_aa_depth_sum, proteome_aa_covered):.6f}",
                "avg_aa_identity": f"{avg_aa_identity:.6f}",
                "avg_aa_identity_weighted": f"{avg_aa_identity_weighted:.6f}",
                "hits_missing_in_gff": hits_missing_here,
            }
        )

    with open(output_tsv, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_out)

    click.echo(f"Wrote: {output_tsv}")
    click.echo(
        f"Blast rows retained: {blast_rows_kept}/{blast_rows_total} "
        f"(unmapped subject proteins: {blast_rows_unmapped_subject})",
        err=True,
    )
    click.echo(
        f"Genomes summarized: {len(rows_out)} "
        f"(missing GFF3: {missing_gff}, hit proteins missing in GFF3: {missing_protein_in_gff}, "
        f"filtered by virus protein coverage: {filtered_by_protein_coverage}, "
        f"filtered by weighted identity: {filtered_by_weighted_identity})",
        err=True,
    )


if __name__ == "__main__":
    cli()