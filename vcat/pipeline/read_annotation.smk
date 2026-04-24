import os
import re

configfile: "config.yaml"

# =============================================================================
# Configuration & Globals
# =============================================================================

DBDIR  = config["database_dir"]
OUTDIR = config["output_dir"]
TMPDIR = config.get("tmpdir", f"{OUTDIR}/tmp")

# Database / reference paths
REF_GENOMES      = f"{DBDIR}/VMR_latest/genomes.fna"
REF_GFF_PATH     = f"{DBDIR}/VMR_latest/gff/"
PROTEIN_2_GENOME = f"{DBDIR}/VMR_latest/genome2protein"


# -----------------------------------------------------------------------------
# Multi-sample configuration
# Accepts:
#   - dict:  {"sample_id": "in=/path/R1.fastq.gz,in2=/path/R2.fastq.gz", ...}
#   - list:  ["in=/path/R1.fastq.gz,in2=/path/R2.fastq.gz", ...]
#   - str:   "in=/path/R1.fastq.gz"  (or with in2=)
# -----------------------------------------------------------------------------
def parse_infiles(infiles_str):
    """Parse 'in=<R1>[,in2=<R2>]' strings into (r1, r2_or_None)."""
    r1 = r2 = None
    for part in str(infiles_str).split(","):
        part = part.strip()
        if part.startswith("in="):
            r1 = part[len("in="):]
        elif part.startswith("in1="):
            r1 = part[len("in1="):]
        elif part.startswith("in2="):
            r2 = part[len("in2="):]
    if not r1:
        raise ValueError("Each sample entry must contain in=<R1> (or in1=<R1>).")
    return r1, r2


def sample_from_fastq(path):
    """Derive a sample name from a FASTQ filename."""
    name = os.path.basename(path)
    name = re.sub(r"\.(fastq|fq)(\.gz)?$", "", name)   # drop .fastq/.fq(.gz)
    name = re.sub(r"([._-])R?[12]$", "", name)         # drop _R1/_R2/.1/.2 etc
    return name


samples_raw = config["sample"]

if isinstance(samples_raw, dict):
    SAMPLES = {k: parse_infiles(v) for k, v in samples_raw.items()}

elif isinstance(samples_raw, list):
    SAMPLES = {}
    for entry in samples_raw:
        r1, r2 = parse_infiles(entry)
        name = sample_from_fastq(r1)
        if name in SAMPLES:
            raise ValueError(
                f"Duplicate sample name '{name}' derived from input list. "
                f"Use a dict in config to assign explicit unique names."
            )
        SAMPLES[name] = (r1, r2)

else:
    r1, r2 = parse_infiles(samples_raw)
    name = sample_from_fastq(r1)
    SAMPLES = {name: (r1, r2)}

SAMPLE_NAMES = list(SAMPLES.keys())


# =============================================================================
# Helper functions for dynamic inputs / params
# =============================================================================

def get_r1(wildcards):
    return SAMPLES[wildcards.sample][0]


def get_r2(wildcards):
    return SAMPLES[wildcards.sample][1]


def get_reads(wildcards):
    """Return a list of existing read files for dependency tracking."""
    r1, r2 = SAMPLES[wildcards.sample]
    return [r1, r2] if r2 else [r1]


def bbmap_in_params(wildcards):
    """Build BBMap-style input arguments (SE vs PE)."""
    r1, r2 = SAMPLES[wildcards.sample]
    if r2:
        return f"in1={r1} in2={r2}"
    return f"in={r1}"


def tool_args(key):
    """Safely fetch optional tool arguments from config."""
    return str(config.get(key) or "").strip()


# =============================================================================
# Target Rule
# =============================================================================

rule all:
    input:
        bam      = expand(f"{OUTDIR}/{{sample}}.sorted.bam",      sample=SAMPLE_NAMES),
        bai      = expand(f"{OUTDIR}/{{sample}}.sorted.bam.bai",  sample=SAMPLE_NAMES),
        pileup   = expand(f"{OUTDIR}/{{sample}}.pileup.tsv",      sample=SAMPLE_NAMES),
        nuc      = expand(f"{OUTDIR}/{{sample}}.vcat_nuc.tsv",    sample=SAMPLE_NAMES),
        diamond  = expand(f"{OUTDIR}/{{sample}}.diamond.m8",       sample=SAMPLE_NAMES),
        protein  = expand(f"{OUTDIR}/{{sample}}.vcat_protein.tsv", sample=SAMPLE_NAMES),


# =============================================================================
# Nucleotide mapping branch
# =============================================================================

rule map_reads:
    input:
        ref   = REF_GENOMES,
        reads = get_reads,
    output:
        sam = temp(f"{TMPDIR}/{{sample}}.sam")
    params:
        bbmap_in       = bbmap_in_params,
        bbmap_index_path = TMPDIR,
        bbmap_args     = lambda wc: config.get("bbmap_args", ""),
    threads: max(1, int(workflow.cores * 0.75))
    log:
        f"{OUTDIR}/logs/{{sample}}_map_reads.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}_map_reads.tsv"
    shell:
        r"""
        mkdir -p {params.bbmap_index_path}
        bbmap.sh ref={input.ref} {params.bbmap_in} out={output.sam} path={params.bbmap_index_path} {params.bbmap_args} &> {log}
        """


rule sort_bam:
    input:
        sam = f"{TMPDIR}/{{sample}}.sam"
    output:
        bam = f"{OUTDIR}/{{sample}}.sorted.bam"
    threads: max(1, int(workflow.cores * 0.75))
    log:
        f"{OUTDIR}/logs/{{sample}}_sort_bam.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}_sort_bam.tsv"
    shell:
        r"""
        samtools sort -@ {threads} -o {output.bam} {input.sam} &> {log}
        """


rule index_bam:
    input:
        bam = f"{OUTDIR}/{{sample}}.sorted.bam"
    output:
        bai = f"{OUTDIR}/{{sample}}.sorted.bam.bai"
    log:
        f"{OUTDIR}/logs/{{sample}}_index_bam.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}_index_bam.tsv"
    shell:
        r"""
        samtools index {input.bam} {output.bai} &> {log}
        """


rule pileup:
    input:
        bam = f"{OUTDIR}/{{sample}}.sorted.bam"
    output:
        tsv = f"{OUTDIR}/{{sample}}.pileup.tsv"
    params:
        pileup_args = lambda wc: tool_args("pileup_args")
    log:
        f"{OUTDIR}/logs/{{sample}}_pileup.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}_pileup.tsv"
    shell:
        r"""
        pileup.sh in={input.bam} out={output.tsv} {params.pileup_args} &> {log}
        """


rule summarize_nuc_mapping:
    input:
        tsv = f"{OUTDIR}/{{sample}}.pileup.tsv"
    output:
        tsv = f"{OUTDIR}/{{sample}}.vcat_nuc.tsv"
    params:
        database_dir = DBDIR,
        summary_args = lambda wc: tool_args("summary_args")
    log:
        f"{OUTDIR}/logs/{{sample}}_summarize_nuc.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}_summarize_nuc.tsv"
    shell:
        r"""
        read_mapping_summary.py -i {input.tsv} -o {output.tsv} -d {params.database_dir} {params.summary_args} &> {log}
        """


# =============================================================================
# Protein mapping branch
# =============================================================================

rule reads_to_fasta:
    input:
        r1 = get_r1,
    params:
        r2 = lambda wc: SAMPLES[wc.sample][1] or "",
    output:
        temp(f"{TMPDIR}/{{sample}}.fasta")
    log:
        f"{OUTDIR}/logs/{{sample}}_reads_to_fasta.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}_reads_to_fasta.tsv"
    shell:
        r"""
        if [ -n "{params.r2}" ]; then
            seqkit fq2fa {input.r1} {params.r2} > {output} 2> {log}
        else
            seqkit fq2fa {input.r1} > {output} 2> {log}
        fi
        """


rule map_reads_to_proteins:
    input:
        fasta = f"{TMPDIR}/{{sample}}.fasta"
    output:
        m8 = f"{OUTDIR}/{{sample}}.diamond.m8"
    params:
        db           = f"{DBDIR}/VMR_latest/diamond_proteins/diamond_proteins.dmnd",
        diamond_args = lambda wc: config.get("diamond_args", ""),
    threads: max(1, int(workflow.cores * 0.75))
    log:
        f"{OUTDIR}/logs/{{sample}}_diamond.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}_diamond.tsv"
    shell:
        r"""
        diamond blastx --query {input.fasta} --out {output.m8} --threads {threads} --db {params.db} \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids slineages \
            {params.diamond_args} &> {log}
        """


rule summarize_protein_mapping:
    input:
        m8 = f"{OUTDIR}/{{sample}}.diamond.m8"
    output:
        tsv = f"{OUTDIR}/{{sample}}.vcat_protein.tsv"
    params:
        ref_gff_path          = REF_GFF_PATH,
        genome2protein          = PROTEIN_2_GENOME,
        min_avg_identity      = config.get("min_avg_aa_identity_weighted", "0.5"),
        min_coverage          = config.get("min_virus_protein_coverage", "0.2"),
        protein_attr          = config.get("protein_attr", "genome_id,protein_id"),
    threads: 1
    log:
        f"{OUTDIR}/logs/{{sample}}_summarize_protein.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}_summarize_protein.tsv"
    shell:
        r"""
        protein_read_mappig_summary.py \
            --blast-columns qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,staxids,slineages \
            --min-avg-aa-identity-weighted {params.min_avg_identity} \
            --min-virus-protein-coverage {params.min_coverage} \
            --protein-attr {params.protein_attr} \
            {input.m8} {params.genome2protein} {params.ref_gff_path} {output.tsv} &> {log}
        """