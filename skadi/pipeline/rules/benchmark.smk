from pathlib import Path
import multiprocessing
import subprocess
import os
from glob import glob
from pathlib import Path
# Extract global variables from the config
DBDIR = config["database_dir"]
OUTDIR = config["results"]
API_PARAMS = config["api"]
AAI_PARAMS = config["aai"]
ANI_PARAMS = config["ani"]
LEVEOUT_LEVEL = config["level"]
BATCH = config["batch"]

# Define database paths
PROTDB = f"{DBDIR}/VMR_latest/mmseqs_proteins/mmseqs_proteins"
GENOMEDB = f"{DBDIR}/VMR_latest/mmseqs_genomes/mmseqs_genomes"
PROFDB = f"{DBDIR}/VMR_latest/mmseqs_pprofiles/mmseqs_pprofiles"
SAMPLD = Path(glob(f"{OUTDIR}/nuc/*_genome.m8")[0]).name.rstrip("_genome.m8")

# Rule to define final output
rule all:
    input:
        expand(f"{OUTDIR}/nuc/{{sample}}_leaveout_{LEVEOUT_LEVEL}_ani.tsv",
        sample=SAMPLD),
        expand(f"{OUTDIR}/prot/{{sample}}_leaveout_{LEVEOUT_LEVEL}_aai.tsv",
        sample=SAMPLD),
        expand(f"{OUTDIR}/prof/{{sample}}_leaveout_{LEVEOUT_LEVEL}_api.tsv",
        sample=SAMPLD),
        expand(f"{OUTDIR}/results/{{sample}}_leaveout_{LEVEOUT_LEVEL}.tsv",
        sample=SAMPLD)

rule cal_ani:
    input:
        f"{OUTDIR}/nuc/{{sample}}_genome.m8"
    output:
        f"{OUTDIR}/nuc/{{sample}}_leaveout_{LEVEOUT_LEVEL}_ani.tsv"
    log:
        f"{OUTDIR}/logs/{{sample}}_leaveout_ani.log"
    params:
        batch = BATCH,
        ani = ANI_PARAMS,
        db = DBDIR,
        level = LEVEOUT_LEVEL
    shell:
        """
        skadi utils ani -i {input} -o {output} \
        --header query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage \
        --level {params.level} --dbdir {params.db} --batch {params.batch} {params.ani} &> {log}
        """

rule cal_aai:
    input:
        m8 = f"{OUTDIR}/prot/{{sample}}_prot.m8",
        gff = f"{OUTDIR}/prot/{{sample}}.gff"
    output:
        f"{OUTDIR}/prot/{{sample}}_leaveout_{LEVEOUT_LEVEL}_aai.tsv"
    log:
        f"{OUTDIR}/logs/{{sample}}_leaveout_aai.log"
    params:
        batch = BATCH,
        aai = AAI_PARAMS,
        db = DBDIR,
        level = LEVEOUT_LEVEL,
    shell:
        """
        skadi utils aai -i {input.m8} -g {input.gff} -d {params.db} -o {output} --topk 300\
        --header query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage \
        --level {params.level} --dbdir {params.db} --batch {params.batch} {params.aai} &> {log}
        """

rule cal_api:
    input:
        m8 = f"{OUTDIR}/prof/{{sample}}_prof.m8",
        gff = f"{OUTDIR}/prot/{{sample}}.gff"
    output:
        f"{OUTDIR}/prof/{{sample}}_leaveout_{LEVEOUT_LEVEL}_api.tsv"
    log:
        f"{OUTDIR}/logs/{{sample}}_leaveout_api.log"
    params:
        batch = BATCH,
        api = API_PARAMS,
        db = DBDIR,
        level = LEVEOUT_LEVEL,
    shell:
        """
        skadi utils api -i {input.m8} -g {input.gff} -d {params.db} -o {output} --topk 300 \
        --header query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage \
        --level {params.level} --dbdir {params.db}  --batch {params.batch} {params.api} &> {log}
        """

rule summarize:
    input:
        DBDIR = DBDIR,
        PROF = f"{OUTDIR}/prof/{{sample}}_leaveout_{LEVEOUT_LEVEL}_api.tsv",
        PROT = f"{OUTDIR}/prot/{{sample}}_leaveout_{LEVEOUT_LEVEL}_aai.tsv",
        NUC = f"{OUTDIR}/nuc/{{sample}}_leaveout_{LEVEOUT_LEVEL}_ani.tsv",

    output:
        f"{OUTDIR}/results/{{sample}}_leaveout_{LEVEOUT_LEVEL}.tsv"
    log:
        f"{OUTDIR}/logs/{{sample}}_leaveout_summarize.log"

    threads: int(workflow.cores * 0.75)

    params:
        ani=ANI_PARAMS,
        api=API_PARAMS,
        aai=AAI_PARAMS

    shell:
        """ 
        postprocess.py {input.DBDIR} {input.NUC} {input.PROT} {input.PROF} {output} {params.ani} {params.aai} {params.api} &> {log}
        
        """