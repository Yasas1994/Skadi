from pathlib import Path
import os

# Extract global variables from the config
DBDIR = config["database_dir"]
OUTDIR = config["results"]
API_PARAMS = config["api"]
AAI_PARAMS = config["aai"]
ANI_PARAMS = config["ani"]
LEVEOUT_LEVEL = config["level"]
BATCH = config["batch"]
METHOD = config.get("method", "cascade")
MIN_CONFIDENCE = config.get("min_confidence", "0.5")

# Define database paths
PROTDB = f"{DBDIR}/VMR_latest/mmseqs_proteins/mmseqs_proteins"
GENOMEDB = f"{DBDIR}/VMR_latest/mmseqs_genomes/mmseqs_genomes"
PROFDB = f"{DBDIR}/VMR_latest/mmseqs_pprofiles/mmseqs_pprofiles"

# Use a checkpoint to discover samples dynamically instead of glob() at parse time.
checkpoint discover_samples:
    output:
        touch(f"{OUTDIR}/.benchmark_samples_discovered")
    params:
        outdir=OUTDIR,
    shell:
        """
        mkdir -p {params.outdir}/nuc
        # If no m8 files exist yet, create a sentinel so the workflow
        # can still proceed; Snakemake will re-evaluate after the
        # checkpoint completes.
        """


def _get_samples(wildcards):
    """Helper to resolve sample names after the checkpoint completes."""
    checkpoint_output = checkpoints.discover_samples.get(**wildcards).output[0]
    nuc_dir = Path(OUTDIR) / "nuc"
    m8_files = list(nuc_dir.glob("*_genome.m8"))
    if not m8_files:
        # Return a dummy sample so Snakemake doesn't crash at DAG build.
        # The actual error will surface when the expected input is missing.
        return ["NO_SAMPLES_FOUND"]
    return [p.name.replace("_genome.m8", "") for p in m8_files]


# Rule to define final output
rule all:
    input:
        lambda wc: expand(
            f"{OUTDIR}/nuc/{{sample}}_leaveout_{LEVEOUT_LEVEL}_ani.tsv",
            sample=_get_samples(wc),
        ),
        lambda wc: expand(
            f"{OUTDIR}/prot/{{sample}}_leaveout_{LEVEOUT_LEVEL}_aai.tsv",
            sample=_get_samples(wc),
        ),
        lambda wc: expand(
            f"{OUTDIR}/prof/{{sample}}_leaveout_{LEVEOUT_LEVEL}_api.tsv",
            sample=_get_samples(wc),
        ),
        lambda wc: expand(
            f"{OUTDIR}/results/{{sample}}_leaveout_{LEVEOUT_LEVEL}.tsv",
            sample=_get_samples(wc),
        ),


rule cal_ani:
    input:
        f"{OUTDIR}/nuc/{{sample}}_genome.m8",
    output:
        f"{OUTDIR}/nuc/{{sample}}_leaveout_{LEVEOUT_LEVEL}_ani.tsv",
    log:
        f"{OUTDIR}/logs/{{sample}}_leaveout_ani.log",
    params:
        batch=BATCH,
        ani=ANI_PARAMS,
        db=DBDIR,
        level=LEVEOUT_LEVEL,
    shell:
        """
        skadi utils ani -i {input} -o {output} \
        --header query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage \
        --level {params.level} --dbdir {params.db} --batch {params.batch} {params.ani} &> {log}
        """


rule cal_aai:
    input:
        m8=f"{OUTDIR}/prot/{{sample}}_prot.m8",
        gff=f"{OUTDIR}/prot/{{sample}}.gff",
    output:
        f"{OUTDIR}/prot/{{sample}}_leaveout_{LEVEOUT_LEVEL}_aai.tsv",
    log:
        f"{OUTDIR}/logs/{{sample}}_leaveout_aai.log",
    params:
        batch=BATCH,
        aai=AAI_PARAMS,
        db=DBDIR,
        level=LEVEOUT_LEVEL,
    shell:
        """
        skadi utils aai -i {input.m8} -g {input.gff} -d {params.db} -o {output} --topk 300 \
        --header query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage \
        --level {params.level} --dbdir {params.db} --batch {params.batch} {params.aai} &> {log}
        """


rule cal_api:
    input:
        m8=f"{OUTDIR}/prof/{{sample}}_prof.m8",
        gff=f"{OUTDIR}/prot/{{sample}}.gff",
    output:
        f"{OUTDIR}/prof/{{sample}}_leaveout_{LEVEOUT_LEVEL}_api.tsv",
    log:
        f"{OUTDIR}/logs/{{sample}}_leaveout_api.log",
    params:
        batch=BATCH,
        api=API_PARAMS,
        db=DBDIR,
        level=LEVEOUT_LEVEL,
    shell:
        """
        skadi utils api -i {input.m8} -g {input.gff} -d {params.db} -o {output} --topk 300 \
        --header query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage \
        --level {params.level} --dbdir {params.db} --batch {params.batch} {params.api} &> {log}
        """


rule summarize:
    input:
        DBDIR=DBDIR,
        PROF=f"{OUTDIR}/prof/{{sample}}_leaveout_{LEVEOUT_LEVEL}_api.tsv",
        PROT=f"{OUTDIR}/prot/{{sample}}_leaveout_{LEVEOUT_LEVEL}_aai.tsv",
        NUC=f"{OUTDIR}/nuc/{{sample}}_leaveout_{LEVEOUT_LEVEL}_ani.tsv",
    output:
        f"{OUTDIR}/results/{{sample}}_leaveout_{LEVEOUT_LEVEL}.tsv",
    log:
        f"{OUTDIR}/logs/{{sample}}_leaveout_summarize.log",
    threads: int(workflow.cores * 0.75)
    params:
        ani=ANI_PARAMS,
        api=API_PARAMS,
        aai=AAI_PARAMS,
        method=METHOD,
        min_confidence=MIN_CONFIDENCE,
    shell:
        """
        postprocess.py {input.DBDIR} {input.NUC} {input.PROT} {input.PROF} {output} \
            --method {params.method} --min-confidence {params.min_confidence} \
            {params.ani} {params.aai} {params.api} &> {log}
        """
