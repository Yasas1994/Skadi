from pathlib import Path
import os

configfile: "config.yaml"

def _detect_system_mem_mb():
    """Return total physical memory in MiB. Platform-aware."""
    try:
        # Best: psutil (add to your env or conda env if available)
        import psutil
        return int(psutil.virtual_memory().total // 1_048_576)
    except Exception:
        pass

    # Linux: /proc/meminfo
    if os.path.exists("/proc/meminfo"):
        try:
            with open("/proc/meminfo") as fh:
                for line in fh:
                    if line.startswith("MemTotal:"):
                        kb = int(line.split()[1])
                        return kb // 1024
        except Exception:
            pass

    # macOS: sysctl
    if shutil.which("sysctl"):
        try:
            out = subprocess.run(
                ["sysctl", "-n", "hw.memsize"],
                capture_output=True, text=True, check=True
            )
            return int(out.stdout.strip()) // 1_048_576
        except Exception:
            pass

    # Windows: wmic
    if shutil.which("wmic"):
        try:
            out = subprocess.run(
                ["wmic", "computersystem", "get", "totalphysicalmemory", "/value"],
                capture_output=True, text=True, check=True
            )
            for line in out.stdout.splitlines():
                if line.strip().startswith("TotalPhysicalMemory="):
                    return int(line.split("=", 1)[1]) // 1_048_576
        except Exception:
            pass

    # Ultimate fallback
    return 8_196


def get_mem_mb(wildcards=None, attempt=1):
    """
    Resolve per-rule memory request (MiB).

    Priority:
      1. User config:   config["mem_mb"]
      2. Environment:   $SKADI_MEM_MB
      3. Auto-detect:   75 % of total system RAM
      4. Hard fallback: 8_196 MiB
    """
    # 1. Explicit user config
    user_mb = config.get("mem_mb")
    if user_mb is not None:
        return int(user_mb)

    # 2. Environment override
    env_mb = os.environ.get("SKADI_MEM_MB")
    if env_mb is not None:
        return int(env_mb)

    # 3. Auto-detect
    total_mb = _detect_system_mem_mb()
    return int(total_mb * 0.75)

# =============================================================================
# Configuration & Globals
# =============================================================================

DBDIR      = config["database_dir"]
OUTDIR     = config["output_dir"]
API_PARAMS = config["api"]
AAI_PARAMS = config["aai"]
ANI_PARAMS = config["ani"]
NUC_SEARCH = config["nuc_search"]

# Database paths
PROTDB   = f"{DBDIR}/VMR_latest/mmseqs_proteins/mmseqs_proteins"
GENOMEDB = f"{DBDIR}/VMR_latest/mmseqs_genomes/mmseqs_genomes"
PROFDB   = f"{DBDIR}/VMR_latest/mmseqs_pprofiles/mmseqs_pprofiles"

# -----------------------------------------------------------------------------
# Multi-sample configuration
# Accepts:
#   - dict: {"sample_id": "/path/to/file.fasta", ...}
#   - list: ["/path/to/a.fasta", "/path/to/b.fa", ...]  (names derived from stems)
#   - str:  "/path/to/single.fasta"
# -----------------------------------------------------------------------------
samples_raw = config["sample"]

if isinstance(samples_raw, dict):
    SAMPLES = samples_raw
elif isinstance(samples_raw, list):
    # Warn if duplicate stems would collide
    seen = set()
    for p in samples_raw:
        stem = Path(p).stem
        if stem in seen:
            raise ValueError(
                f"Duplicate sample name '{stem}' derived from input list. "
                f"Use a dict in config to assign explicit unique names."
            )
        seen.add(stem)
    SAMPLES = {Path(p).stem: p for p in samples_raw}
else:
    SAMPLES = {Path(samples_raw).stem: samples_raw}

SAMPLE_NAMES = list(SAMPLES.keys())


def get_sample_path(wildcards):
    """Resolve the original FASTA path for a given sample wildcard."""
    return SAMPLES[wildcards.sample]


# =============================================================================
# Target Rule
# =============================================================================

rule all:
    input:
        results  = expand(f"{OUTDIR}/results/{{sample}}.tsv", sample=SAMPLE_NAMES),
        nuc_m8   = expand(f"{OUTDIR}/nuc/{{sample}}_genome.m8", sample=SAMPLE_NAMES),
        nuc_ani  = expand(f"{OUTDIR}/nuc/{{sample}}_genome_ani.tsv", sample=SAMPLE_NAMES),
        prot_m8  = expand(f"{OUTDIR}/prot/{{sample}}_prot.m8", sample=SAMPLE_NAMES),
        prot_aai = expand(f"{OUTDIR}/prot/{{sample}}_prot_aai.tsv", sample=SAMPLE_NAMES),
        prof_m8  = expand(f"{OUTDIR}/prof/{{sample}}_prof.m8", sample=SAMPLE_NAMES),
        prof_api = expand(f"{OUTDIR}/prof/{{sample}}_prof_api.tsv", sample=SAMPLE_NAMES),


# =============================================================================
# Nucleotide Annotation (runs on original input)
# =============================================================================

rule annot_nuc:
    input:
        get_sample_path
    output:
        m8  = f"{OUTDIR}/nuc/{{sample}}_genome.m8",
        tmp = temp(directory(f"{OUTDIR}/tmp/nuc/{{sample}}")),
    log:
        f"{OUTDIR}/logs/{{sample}}_annot_nuc.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}_annot_nuc.tsv"
    params:
        db          = GENOMEDB,
        search_type = 3 if NUC_SEARCH == "mmseqs_blastn" else 2,
    threads: max(1, int(workflow.cores * 0.75))
    resources:
        mem_mb = get_mem_mb
    shell:
        """
        mkdir -p {output.tmp}
        mmseqs easy-search {input} {params.db} {output.m8} {output.tmp} \
            --threads {threads} \
            -s 7 \
            --search-type {params.search_type} \
            --format-output query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage \
            &> {log}
        """


# =============================================================================
# Checkpoint: split fasta dynamically per sample
# =============================================================================

checkpoint split_fasta:
    input:
        get_sample_path
    output:
        dir = directory(f"{OUTDIR}/tmp/{{sample}}_parts")
    params:
        max_splits = max(1, workflow.cores),
    threads: 4
    log:
        f"{OUTDIR}/logs/{{sample}}_split_fasta.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}_split_fasta.tsv"
    shell:
        """
        mkdir -p {output.dir}

        # Count sequences at runtime (no parse-time subprocess!)
        NUM_SEQS=$(seqkit stats {input} | awk 'NR==2{{print $4}}' | tr -d ',')
        SPLITS=$NUM_SEQS
        if [ $SPLITS -gt {params.max_splits} ]; then
            SPLITS={params.max_splits}
        fi
        if [ $SPLITS -lt 1 ]; then
            SPLITS=1
        fi

        # Run seqkit into a staging directory
        STAGING="{output.dir}/.staging"
        mkdir -p $STAGING
        seqkit split2 --by-part $SPLITS --threads {threads} {input} -O $STAGING &> {log}

        # Normalize seqkit output to predictable names: part_001.fasta, part_002.fasta, ...
        i=1
        for f in $STAGING/*; do
            [ -e "$f" ] || continue
            num=$(printf "%03d" "$i")
            mv "$f" "{output.dir}/part_${{num}}.fasta"
            i=$((i + 1))
        done
        rm -rf $STAGING
        """


def get_split_parts(wildcards):
    """Discover part indices from checkpoint output for a given sample."""
    checkpoint_dir = checkpoints.split_fasta.get(**wildcards).output.dir
    parts = glob_wildcards(os.path.join(checkpoint_dir, "part_{part}.fasta")).part
    return parts


# =============================================================================
# Protein Prediction
# =============================================================================

rule run_prodigal:
    input:
        f"{OUTDIR}/tmp/{{sample}}_parts/part_{{part}}.fasta"
    output:
        faa = temp(f"{OUTDIR}/tmp/{{sample}}_{{part}}.faa"),
        gff = temp(f"{OUTDIR}/tmp/{{sample}}_{{part}}.gff"),
    log:
        f"{OUTDIR}/logs/{{sample}}_{{part}}_run_prodigal.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}_{{part}}_run_prodigal.tsv"
    threads: 1
    shell:
        """
        prodigal-gv -i {input} -a {output.faa} -o {output.gff} -p meta -f gff &> {log}
        """


rule merge_prodigal:
    input:
        faa = lambda wc: expand(
            f"{OUTDIR}/tmp/{wc.sample}_{{part}}.faa",
            part=get_split_parts(wc),
        ),
        gff = lambda wc: expand(
            f"{OUTDIR}/tmp/{wc.sample}_{{part}}.gff",
            part=get_split_parts(wc),
        ),
    output:
        faa = f"{OUTDIR}/prot/{{sample}}.faa",
        gff = f"{OUTDIR}/prot/{{sample}}.gff",
    log:
        f"{OUTDIR}/logs/{{sample}}_merge_prodigal.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}_merge_prodigal.tsv"
    threads: 1
    shell:
        """
        cat {input.gff} > {output.gff} 2> {log}
        cat {input.faa} > {output.faa} 2> {log}
        """


# =============================================================================
# Protein Annotation
# =============================================================================

rule annot_prot:
    input:
        f"{OUTDIR}/prot/{{sample}}.faa"
    output:
        m8  = f"{OUTDIR}/prot/{{sample}}_prot.m8",
        tmp = temp(directory(f"{OUTDIR}/tmp/prot/{{sample}}")),
    log:
        f"{OUTDIR}/logs/{{sample}}_annot_prot.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}_annot_prot.tsv"
    params:
        db = PROTDB,
    threads: max(1, int(workflow.cores * 0.75))
    resources:
        mem_mb = get_mem_mb
    shell:
        """
        mkdir -p {output.tmp}
        mmseqs easy-search {input} {params.db} {output.m8} {output.tmp} \
            --threads {threads} \
            -s 6 \
            --split-memory-limit {resources.mem_mb}M \
            --format-output query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage \
            &> {log}
        """


# =============================================================================
# Profile Annotation
# =============================================================================

rule annot_prof_search:
    input:
        faa = f"{OUTDIR}/prot/{{sample}}.faa",
        db  = PROFDB,          # e.g. .../mmseqs_pprofiles
    output:
        m8 = f"{OUTDIR}/prof/{{sample}}_prof.m8",
    params:
        tmp_dir      = f"{OUTDIR}/tmp/prof/{{sample}}_search",
        query_db     = f"{OUTDIR}/tmp/prof/{{sample}}_search/query_db",
        result_db    = f"{OUTDIR}/tmp/prof/{{sample}}_search/result_db",
        max_accept   = 100,
        cov          = 0.1,
        start_sens   = 1,
        s            = 7,
        sens_steps   = 3,
        prefilter    = 0, # 2: nofilter
        fmt          = "query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage",
    threads: max(1, int(workflow.cores * 0.75))
    resources:
        mem_mb = get_mem_mb
    log:
        f"{OUTDIR}/logs/{{sample}}_annot_prof_search.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}_annot_prof_search.tsv"
    shell:
        """
        mkdir -p {params.tmp_dir}

        # 1. createdb
        mmseqs createdb {input.faa} {params.query_db} --dbtype 1 >> {log} 2>&1

        # 2. search
        mmseqs search  {params.query_db} {input.db} {params.result_db} {params.tmp_dir} \
            --threads {threads} \
            -c {params.cov} \
            --start-sens {params.start_sens} \
            -s {params.s} \
            --sens-steps {params.sens_steps} \
            --prefilter-mode {params.prefilter} \
            --max-accept {params.max_accept} \
            --split-memory-limit {resources.mem_mb}M \
            >> {log} 2>&1

        # 3. convertalis
        mmseqs convertalis {params.query_db} {input.db} {params.result_db} {output.m8} \
            --threads {threads} \
            --format-output {params.fmt} \
            >> {log} 2>&1

        # 4. tidy up intermediate DBs (optional, mimics easy-search)
        mmseqs rmdb {params.query_db}  >> {log} 2>&1 || true
        mmseqs rmdb {params.result_db} >> {log} 2>&1 || true
        """

rule annot_prof_postprocess:
    input:
        m8 = f"{OUTDIR}/prof/{{sample}}_prof.m8",
    output:
        m8 = f"{OUTDIR}/prof/{{sample}}_prof_corrected.m8",  # or overwrite in-place
    params:
        db_root = DBDIR,
    log:
        f"{OUTDIR}/logs/{{sample}}_annot_prof_postprocess.log"
    shell:
        """
        cp {input.m8} {output.m8}
        correct_profile_taxinfo.py -d {params.db_root} -i {output.m8} &> {log}
        """

# =============================================================================
# Distance Calculations
# =============================================================================

rule cal_ani:
    input:
        f"{OUTDIR}/nuc/{{sample}}_genome.m8"
    output:
        f"{OUTDIR}/nuc/{{sample}}_genome_ani.tsv"
    log:
        f"{OUTDIR}/logs/{{sample}}_ani.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}_ani.tsv"
    params:
        batch = 1000,
        ani   = ANI_PARAMS,
    shell:
        """
        skadi utils ani -i {input} \
            --header query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage \
            --batch {params.batch} {params.ani} \
            > {output} 2> {log}
        """


rule cal_aai:
    input:
        m8  = f"{OUTDIR}/prot/{{sample}}_prot.m8",
        gff = f"{OUTDIR}/prot/{{sample}}.gff"
    output:
        f"{OUTDIR}/prot/{{sample}}_prot_aai.tsv"
    log:
        f"{OUTDIR}/logs/{{sample}}_aai.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}_aai.tsv"
    params:
        batch = 1000,
        aai   = AAI_PARAMS,
        db    = DBDIR,
    shell:
        """
        skadi utils aai -i {input.m8} -g {input.gff} -d {params.db} \
            --header query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage \
            --batch {params.batch} {params.aai} \
            > {output} 2> {log}
        """


rule cal_api:
    input:
        m8  = f"{OUTDIR}/prof/{{sample}}_prof.m8",
        gff = f"{OUTDIR}/prot/{{sample}}.gff"
    output:
        f"{OUTDIR}/prof/{{sample}}_prof_api.tsv"
    log:
        f"{OUTDIR}/logs/{{sample}}_api.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}_api.tsv"
    params:
        batch = 1000,
        api   = API_PARAMS,
        db    = DBDIR,
    shell:
        """
        skadi utils api -i {input.m8} -g {input.gff} -d {params.db} \
            --header query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage \
            --batch {params.batch} {params.api} \
            > {output} 2> {log}
        """


# =============================================================================
# Summarization
# =============================================================================

rule summarize:
    input:
        prof = f"{OUTDIR}/prof/{{sample}}_prof_api.tsv",
        prot = f"{OUTDIR}/prot/{{sample}}_prot_aai.tsv",
        nuc  = f"{OUTDIR}/nuc/{{sample}}_genome_ani.tsv",
    output:
        f"{OUTDIR}/results/{{sample}}.tsv"
    log:
        f"{OUTDIR}/logs/{{sample}}_summarize.log"
    benchmark:
        f"{OUTDIR}/benchmarks/{{sample}}_summarize.tsv"
    threads: max(1, int(workflow.cores * 0.75))
    params:
        db  = DBDIR,
        ani = ANI_PARAMS,
        api = API_PARAMS,
        aai = AAI_PARAMS,
    shell:
        """
        postprocess.py {params.db} {input.nuc} {input.prot} {input.prof} {output} {params.ani} {params.aai} {params.api} &> {log}
        """