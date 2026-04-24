import click
import subprocess
import multiprocessing
import os
import sys
import yaml
import polars as pl
from pathlib import Path
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
from vcat.utils import ani_summary, axi_summary, index_m8, load_chunk
from .color_logger import logger
from importlib.metadata import version, PackageNotFoundError


# Define the directory containing the pipeline files
PIPELINE_DIR = os.path.join(os.path.dirname(__file__), "./pipeline")
CONFIG = os.path.join(PIPELINE_DIR, "config.yaml")
HEADER = "query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage"
CONFIG_CONTENT = yaml.safe_load(open(CONFIG, "r"))

__version__ = ""
try:
    __version__ = version("vcat")
except PackageNotFoundError:
    # package is not installed
    pass
def load_configfile(file_path):
    with open(file_path, "r") as f:
        return yaml.safe_load(f)

def parse_csv(ctx, param, value):
    """Split the comma-separated input into a list."""
    if value:
        return value.split(",")
    return []

def format_databases(config):
    default_url = ""
    default = ""
    choices = []
    for i in config["downloads"]:
        choices.append(i.get("name"))
        if i.get("latest"):
            default += i.get("name")
            default_url += i.get("link")
    return default, choices, default_url

def handle_max_mem(max_mem, profile):
    "Specify maximum virtual memory to use by atlas."
    "For numbers >1 its the memory in GB. "
    "For numbers <1 it's the fraction of available memory."

    if profile is not None:
        if max_mem is not None:
            logger.info(
                "Memory requirements are handled by the profile, I ignore max-mem argument."
            )
        # memory is handled via the profile, user should know what he is doing
        return ""
    else:
        import psutil
        from math import floor

        # calulate max  system meory in GB (float!)
        max_system_memory = psutil.virtual_memory().total / (1024**3)

        if max_mem is None:
            max_mem = 0.95
        if max_mem > 1:
            if max_mem > max_system_memory:
                logger.critical(
                    f"You specified {max_mem} GB as maximum memory, but your system only has {floor(max_system_memory)} GB"
                )
                sys.exit(1)

        else:
            max_mem = max_mem * max_system_memory

        # specify max_mem_string including java mem and max mem

        return f" --resources mem={floor(max_mem)} mem_mb={floor(max_mem*1024)} java_mem={floor(0.85* max_mem)} "

def get_snakefile(file=f"{PIPELINE_DIR }/Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf

def update_config(config_path, data: dict):
    import yaml

    # Step 1: Read the YAML file
    with open(config_path, "r") as file:
        config = yaml.safe_load(file)  # Load the YAML into a Python dictionary

    # Step 2: Update the configuration
    for k, v in data.items():
        config[k] = v

    # Step 3: Write back to the YAML file
    with open(config_path, "w") as file:
        yaml.safe_dump(config, file, default_flow_style=False)


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context
def cli(obj):
    r"""
    Vcat: a command-line tool-kit for adding ICTV taxonomy annotations to virus contigs,
    mapping reads to virus genomes and much more.
    (https://github.com/Yasas1994/vcat)"""
    pass


@cli.command(
    context_settings=dict(ignore_unknown_options=True, show_default=True),
    short_help="run contig annotation workflow",
    help="""
    The Virus contig annotation tool (Vcat) is a straightforward, homology-based application designed to 
    provide taxonomic annotations for viral contigs. 

    """,
)
@click.option(
    "-i",
    "--input",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="input file/s to run vcat on",
    required=True,
)
@click.option(
    "-o",
    "--output",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="dir to store vcat results",
    required=True,
)
@click.option(
    "-d",
    "--database",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="dir to vcat database",
    required=False,
)
@click.option(
    "--nuc_search",
    type=click.Choice(["blastn", "tbalstx", "mmseqs_blastn", "mmseqs_tblastx"]),
    default="mmseqs_blastn",
    show_default=True,
    help="nucelotide search algorithm (blastn, tblastx, mmseqs_blastn, mmseqs_tblastx).",
)
@click.option(
    "--prot_search",
    type=click.Choice(["blast", "diamond", "mmseqs"]),
    default="mmseqs",
    show_default=True,
    help="protein search algorithm (blast, diamond, mmseqs).",
)
@click.option(
    "--prof_search",
    type=click.Choice(["mmseqs", "hmmer"]),
    default="mmseqs",
    show_default=True,
    help="nucelotide seach algorithm (mmseqs, hmmer).",
)
@click.option(
    "-j",
    "--jobs",
    type=int,
    default=multiprocessing.cpu_count(),
    show_default=True,
    help="use at most this many jobs in parallel (see cluster submission for more details).",
)
@click.option(
    "--tapif",
    type=float,
    default=0.3,
    help="assign sequences above this tapi threshold to families",
    required=False,
)
@click.option(
    "--tapio",
    type=float,
    default=0.15,
    help="assign sequences above this tapi threshold to orders",
    required=False,
)
@click.option(
    "--tapic",
    type=float,
    default=0.15,
    help="assign sequences above this tapi threshold to classes",
    required=False,
)
@click.option(
    "--tapip",
    type=float,
    default=0.15,
    help="assign sequences above this tapi threshold to phyla",
    required=False,
)
@click.option(
    "--tapik",
    type=float,
    default=0.15,
    help="assign sequences above this tapi threshold to kingdoms",
    required=False,
)
@click.option(
    "--taaig",
    type=float,
    default=0.49,
    help="assign sequences above this taai threshold to genera",
    required=False,
)
@click.option(
    "--taaif",
    type=float,
    default=0.3,
    help="assign sequences above this taai threshold to families",
    required=False,
)
@click.option(
    "--taaio",
    type=float,
    default=0.3,
    help="assign sequences above this taai threshold to orders",
    required=False,
)
@click.option(
    "--taaic",
    type=float,
    default=0.3,
    help="assign sequences above this taai threshold to classes",
    required=False,
)
@click.option(
    "--taaip",
    type=float,
    default=0.3,
    help="assign sequences above this taai threshold to phyla",
    required=False,
)
@click.option(
    "--taaik",
    type=float,
    default=0.3,
    help="assign sequences above this taai threshold to kingdoms",
    required=False,
)
@click.option(
    "--tanis",
    type=float,
    default=0.81,
    help="assign sequences above this taai threshold to species",
    required=False,
)
@click.option(
    "--tanig",
    type=float,
    default=0.49,
    help="assign sequences above this taai threshold to genera",
    required=False,
)
@click.option(
    "--batch",
    type=int,
    default=5000,
    help="number of records to process at a time",
    required=False,
)
@click.option(
    "-n",
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="Test execution.",
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def contigs(input, output, database, jobs, batch, snakemake_args, **kwargs):
    """
    Runs Vcat pipeline on contigs

    Most snakemake arguments can be appended to the command for more info see 'snakemake --help'
    """

    logger.info(f"vcat version: {__version__}")
    conf = load_configfile(CONFIG)
    db_dir = database or conf["database_dir"]
    logger.info(f"database version: {Path(db_dir).name}")

    taai_parts, tapi_parts, tani_parts = [], [], []

    for k, v in kwargs.items():
        if v is None:
            continue
        if k.startswith("taai"):
            taai_parts.append(f"--{k} {v}")
        elif k.startswith("tapi"):
            tapi_parts.append(f"--{k} {v}")
        elif k.startswith("tani"):
            tani_parts.append(f"--{k} {v}")

    taai_params = " ".join(taai_parts)
    tapi_params = " ".join(tapi_parts)
    tani_params = " ".join(tani_parts)

    # nuc_search: only pass if explicitly provided
    nuc_search = kwargs.get("nuc_search", None)
    if isinstance(nuc_search, bool):
        nuc_search = str(nuc_search).lower()  # "true"/"false" plays nicer in configs

    cmd_parts = [
        "snakemake",
        f"--snakefile {get_snakefile('./pipeline/contig_annotation.smk')}",
        f"--jobs {jobs}",
        "--rerun-incomplete",
        f"--configfile {CONFIG}",
        "--scheduler greedy",
        "--show-failed-logs",
        "--groups group1=1",
        "--config",
        f"database_dir='{db_dir}'",
        f"sample='{input}'",
        f"output_dir='{output}'",
        f"api='{tapi_params}'",
        f"aai='{taai_params}'",
        f"ani='{tani_params}'",
        f"batch='{batch}'",
    ]

    if nuc_search is not None:
        cmd_parts.append(f"nuc_search='{nuc_search}'")

    if snakemake_args:
        cmd_parts.extend(snakemake_args)

    cmd = " ".join(cmd_parts)

    logger.info(cmd)
    logger.debug("Executing: %s" % cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logger.critical(e)
        exit(1)


@cli.command(
    context_settings=dict(ignore_unknown_options=True, show_default=True),
    short_help="run read annotation workflow",
    help="""
    The Virus contig annotation tool (Vcat) is a straightforward, homology-based application designed to 
    provide taxonomy annotations to virus contigs and mapping reads directly to virus genomes.

    usage (paired-end): vcat reads [OPTIONS] -i1 pair1.fastq -i2 pair2.fastq -o mapping_results.tsv

    usage (single-end): vcat reads [OPTIONS] -i1 pair1.fastq -o mapping_results.tsv
    """,
)
@click.option(
    "-in",
    "--input",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="input read file1 to run vcat on",
    required=True,
)
@click.option(
    "-in2",
    "--input2",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="input read file2 (paired-end) to run vcat on",
    required=False,
)
@click.option(
    "-o",
    "--output",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="dir to store vcat results",
    required=True,
)
@click.option(
    "-d",
    "--database",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="dir to vcat database",
    required=False,
)
@click.option(
    "-j",
    "--jobs",
    type=int,
    default=multiprocessing.cpu_count(),
    show_default=True,
    help="use at most this many jobs in parallel (see cluster submission for more details).",
)
@click.option(
    "--bbmap_args",
    default="nodisk=t minid=0.95 maxindel=2 -Xmx20g",
    help="Extra arguments passed directly to BBMap (e.g. 'minid=0.95 maxindel=3')",
)
@click.option(
    "--pileup-args",
    default=" qtrim=t trimq=10 border=5 secondary=f delcoverage=f minmapq=20",
    help=(
        "Extra arguments passed to pileup.sh / CoveragePileup "
        "(e.g. 'minmapq=20 minbaseq=20 mincov=2 secondary=f delcov=f')"
    ),
)
@click.option(
    "--summary-args",
    default=" -cp 0.5 -af 1 -mtr 100",
    help=(
        "coverage percentage = cp"
        "Average fold = af"
        "Min total reads = mtr"
        "-cp 50.0 af 1 mtr 100"
    ),
)
@click.option(
    "--profile",
    default=None,
    help="snakemake profile e.g. for cluster execution.",
)
@click.option(
    "-n",
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="Test execution.",
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def reads(input, input2, output, database, jobs, profile, dryrun, bbmap_args, pileup_args, summary_args, snakemake_args):
    """
    Runs Vcat pipeline on reads

    Most snakemake arguments can be appended to the command for more info see 'snakemake --help'
    """

    logger.info(f"vcat version: {__version__}")
    conf = load_configfile(CONFIG)
    db_dir = database or conf["database_dir"]
    logger.info(f"database version: {Path(db_dir).name}")

    cmd = [
        "snakemake",
        f"--snakefile {get_snakefile('./pipeline/read_annotation.smk')}",
        f"--jobs {jobs}",
        "--rerun-incomplete",
        f"--configfile {CONFIG}",
        "--scheduler greedy",
        "--show-failed-logs",
        "--groups group1=1",
        f"--config database_dir='{db_dir}' sample='in={input}, in2={input2}' output_dir='{output}' bbmap_args='{bbmap_args}' pileup_args='{pileup_args}' summary_args='{summary_args}'",
    ]

    if dryrun:
        cmd.append("--dry-run")

    if profile:
        cmd.append(f"--profile {profile}")

    if snakemake_args:
        if snakemake_args[0].startswith("-"):
            cmd.extend(snakemake_args)
        else:
            cmd.append("--")
            cmd.extend(snakemake_args)

    cmd = " ".join(cmd)

    logger.debug("Executing: %s" % cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logger.critical(e)
        exit(1)


# Download and build
@cli.command(
    context_settings=dict(ignore_unknown_options=True, show_default=True),
    short_help="download and build reference databases",
)
@click.option(
    "-d",
    "--db-dir",
    help="location to store databases",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    required=True,
)
@click.option(
    "-j",
    "--jobs",
    default=1,
    type=int,
    show_default=True,
    help="number of simultaneous downloads",
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def preparedb(db_dir, jobs, snakemake_args):
    """Executes a snakemake workflow to downlod and building the databases"""
    logger.info("Building taxdb")
    cmd_parts = [
        "snakemake",
        f"--snakefile {get_snakefile('./pipeline/rules/taxdump.smk')}",
        f"--jobs {jobs}",
        "--rerun-incomplete",
        f"--configfile {CONFIG}",
        "--scheduler greedy",
        "--show-failed-logs",
        "--config",
        f"database_dir='{db_dir}'",
    ]

    if snakemake_args:
        if snakemake_args[0].startswith("-"):
            cmd_parts.extend(snakemake_args)
        else:
            cmd_parts.append("--")
            cmd_parts.extend(snakemake_args)

    cmd = " ".join(cmd_parts)

    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logger.critical(e)
        exit(1)

    logger.info("Building mmseqs databases")

    cmd_parts = [
        "snakemake",
        f"--snakefile {get_snakefile('./pipeline/rules/createdb.smk')}",
        f"--jobs {jobs}",
        "--rerun-incomplete",
        f"--configfile {CONFIG}",
        "--scheduler greedy",
        "--show-failed-logs",
        "--config",
        f"database_dir='{db_dir}'",
    ]

    if snakemake_args:
        if snakemake_args[0].startswith("-"):
            cmd_parts.extend(snakemake_args)
        else:
            cmd_parts.append("--")
            cmd_parts.extend(snakemake_args)

    cmd = " ".join(cmd_parts)


    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logger.critical(e)
        exit(1)

    logger.info(f"Adding {db_dir} to config.yaml")
    update_config(config_path=CONFIG, data={"database_dir": db_dir})


# Download pre-built from server
default_db, choices, default_url = format_databases(config=CONFIG_CONTENT)


@cli.command(
    context_settings=dict(ignore_unknown_options=True, show_default=True),
    short_help="pull pre-built databases from a remote server",
)
@click.option(
    "-d",
    "--db-dir",
    help="location to store databases",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    required=True,
)
@click.option(
    "--dbversion",
    default=default_db,
    type=click.Choice(choices, case_sensitive=False),
    show_default=True,
    help="version of the database to download",
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def downloaddb(db_dir, dbversion, snakemake_args):
    "pull pre-built databases from a remote server"

    logger.info(f"downloading {dbversion} database from remote server")
    cmd_parts = [
        "snakemake",
        f"--snakefile {get_snakefile('./pipeline/rules/download.smk')}",
        "--jobs 1",
        "--rerun-incomplete",
        f"--configfile {CONFIG}",
        "--scheduler greedy",
        "--show-failed-logs",
        "--config",
        f"database_dir='{db_dir}'",
        f"dbversion='{dbversion}'",
        f"dburl='{default_url}'",
    ]

    if snakemake_args:
        if snakemake_args[0].startswith("-"):
            cmd_parts.extend(snakemake_args)
        else:
            cmd_parts.append("--")
            cmd_parts.extend(snakemake_args)

    cmd = " ".join(cmd_parts)

    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logger.critical(e)
        exit(1)

    logger.info(f"Adding {db_dir} to config.yaml")
    update_config(
        config_path=CONFIG, data={"database_dir": str(Path(db_dir) / dbversion)}
    )


# utility functions
@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context
def utils(obj):
    """
    tool chain for calculating ani, aai and visualizations
    """
    pass


@utils.command(
    context_settings=dict(ignore_unknown_options=True, show_default=True),
    help="""calculates ani from mmseqs ICTV genome comparision results
            and writes the results to <input>_ani.tsv

            output includes qcov (i.e synonymous to alignment fraction),
            ani (average nucleotide identity) and tani (ani x qcov)

            usage                                                                      
            -----                                                                      

            vcat utils ani [OPTIONS] -i contigs.fasta

        """,
)
@click.option(
    "-i",
    "--input",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="BLAST tabular (m8) output file",
    required=True,
)
@click.option(
    "-o",
    "--output",
    type=str,
    default=None,
    help="output file path",
    required=False,
)
@click.option(
    "--header",
    callback=parse_csv,
    help="columnames of the m8 file",
    default=HEADER,
    required=False,
)
@click.option(
    "--tanig",
    type=float,
    default=0.49,
    help="assign sequences above this taai threshold to genera",
    required=False,
)
@click.option(
    "--tanis",
    type=float,
    default=0.81,
    help="assign sequences above this taai threshold to species",
    required=False,
)
@click.option(
    "--all",
    is_flag=True,
    default=False,
    help="get ani for all hits per query sequence. by default only outputs the besthits",
    required=False,
)
@click.option(
    "--level",
    type=str,
    default=None,
    help="leave out taxonomic level (species, genus, family, class)",
    required=False,
)
@click.option(
    "-d",
    "--dbdir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="path to the ref databases",
    required=False,
)
@click.option(
    "--batch",
    type=int,
    default=5000,
    help="number of records to process at a time",
    required=False,
)
def ani(input, output, header, level, dbdir, all, batch, **kwargs):
    """
    calculates the average nucleotide identity and coverage of a query sequence
    to the best
    """
    if (level or dbdir) and not (level and dbdir):
            logger.error(f"{level} and {dbdir} are mutually inclusive")
            exit(1)


    CHUNK_SIZE = batch
    file_name = os.path.basename(input)

    index = index_m8(input, kind="ani")
    tmp_files = []
    with logging_redirect_tqdm():
        for i in tqdm(range(0, len(index), CHUNK_SIZE), ncols=70, ascii=" ="):
            finput = load_chunk(
                input, index=index, recstart=i, recend=min(i + CHUNK_SIZE, len(index))
            )
            outfile = os.path.join(
                    os.path.dirname(input),
                    f"{os.path.splitext(file_name)[0]}_ani_{min(i + CHUNK_SIZE, len(index))}.tsv",
                )
            status = ani_summary(finput, all=all, header=header, level=level, dbdir=dbdir)
            if isinstance(status, pl.DataFrame):
                status.write_csv(outfile, separator="\t")
                logger.info(f"{outfile} updated")
                tmp_files.append(outfile)

            else:
                # exit code 1
                logger.error("error occured!")
                logger.exception(status)
                exit(1)

    logger.info("trying to merge temporary files")
    tmp = [pl.read_csv(f, separator="\t") for f in tmp_files]
    tmp = [i for i in tmp if not i.is_empty()]
    if not output:
        outfile = os.path.join(
                os.path.dirname(input), f"{os.path.splitext(file_name)[0]}_ani.tsv"
            )
    else:
        outfile = output
    if tmp:
        df = pl.concat(tmp)

        df.write_csv(outfile, separator="\t")
        logger.info(f"{outfile} updated")
    else:
        logger.info("all tables are empty")
        Path(outfile).touch()

    # Remove temporary TSV files
    for file in tmp_files:
        os.remove(file)


@utils.command(
    context_settings=dict(ignore_unknown_options=True, show_default=True),
    help="""
            calculates aai from mmseqs ICTV viral protein comparision results
            and writes the results to <input>_aai.tsv

            output includes qcov (i.e number of genes shared between the query
            and target over number of genes on the query),
            aai (average aminoacid identity) and taai (aai x qcov)

            usage                                                                      
            -----                                                                      

            vcat utils aai  [OPTIONS] -i contigs.fasta -g configs.gff -d [DBDIR]

        """,
)
@click.option(
    "-i",
    "--input",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="BLAST tabular (m8) output file",
    required=True,
)
@click.option(
    "-o",
    "--output",
    type=str,
    default=None,
    help="output file path",
    required=False,
)
@click.option(
    "-g",
    "--gff",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="gff file with gene cordinates of the query sequences",
    required=True,
)
@click.option(
    "--header",
    callback=parse_csv,
    help="columnames of the m8 file ",
    default=HEADER,
    required=False,
)
@click.option(
    "--taaig",
    type=float,
    default=0.3,
    help="assign sequences above this taai threshold to genera",
    required=False,
)
@click.option(
    "--taaif",
    type=float,
    default=0.3,
    help="assign sequences above this taai threshold to families",
    required=False,
)
@click.option(
    "--taaio",
    type=float,
    default=0.3,
    help="assign sequences above this taai threshold to orders",
    required=False,
)
@click.option(
    "--taaic",
    type=float,
    default=0.3,
    help="assign sequences above this taai threshold to classes",
    required=False,
)
@click.option(
    "--taaip",
    type=float,
    default=0.3,
    help="assign sequences above this taai threshold to phyla",
    required=False,
)
@click.option(
    "--taaik",
    type=float,
    default=0.3,
    help="assign sequences above this taai threshold to kingdoms",
    required=False,
)
@click.option(
    "-d",
    "--dbdir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="path to the ref databases",
    required=True,
)
@click.option(
    "--batch",
    type=int,
    default=5000,
    help="number of records to process at a time",
    required=False,
)
@click.option(
    "--topk",
    type=int,
    default=5,
    help="number of hits to select per query",
    required=False,
)
@click.option(
    "--level",
    type=str,
    default=None,
    help="leave out taxonomic level (species, genus, family, class)",
    required=False,
)
@click.option(
    "--all",
    is_flag=True,
    default=False,
    help="get aai for all top-k hits per query sequence. by default only outputs the besthit",
    required=False,
)
def aai(
    input, output, header, taaig, taaif, taaio, taaic, taaip, taaik, batch, dbdir, gff, topk, level, all
):
    """
    calculates the average aminoacid identity and coverage of a query sequence
    to the genomes in the target database
    """
    CHUNK_SIZE = batch
    THRESHOLDS = {
        "genus": taaig,
        "family": taaif,
        "order": taaio,
        "class": taaic,
        "phylum": taaip,
        "kingdom": taaik,
    }
    file_name = os.path.basename(input)

    index = index_m8(input, kind="axi")
    tmp_files = []
    with logging_redirect_tqdm():
        for i in tqdm(range(0, len(index), CHUNK_SIZE), ncols=70, ascii=" ="):
            finput = load_chunk(
                input, index=index, recstart=i, recend=min(i + CHUNK_SIZE, len(index))
            )
            outfile = os.path.join(
                os.path.dirname(input),
                f"{os.path.splitext(file_name)[0]}_aai_{min(i + CHUNK_SIZE, len(index))}.tsv",
            )
            status = axi_summary(
                finput, gff, dbdir, header, THRESHOLDS, top_k=topk, kind="aai", all=all, level=level
            )
            if isinstance(status, pl.DataFrame):
                status.write_csv(outfile, separator="\t")
                logger.info(f"{outfile} updated")
                tmp_files.append(outfile)

            else:
                # exit code 1
                logger.error("error occured!")
                logger.exception(status)
                exit(1)

    logger.info("merging temporary files")
    tmp = [pl.read_csv(f, separator="\t") for f in tmp_files]
    tmp = [i for i in tmp if not i.is_empty()]
    if not output:
        outfile = os.path.join(
                os.path.dirname(input), f"{os.path.splitext(file_name)[0]}_aai.tsv"
            )
    else:
        outfile = output
    if tmp:
        df = pl.concat(tmp)

        df.write_csv(outfile, separator="\t")
        logger.info(f"{outfile} updated")
    else:
        logger.info("all tables are empty")
        Path(outfile).touch()

    # Remove temporary TSV files
    for file in tmp_files:
        os.remove(file)


@utils.command(
    context_settings=dict(ignore_unknown_options=True, show_default=True,),
    help="""
            calculates api from mmseqs ICTV viral protein comparision results
            and writes the results to <input>_api.tsv

            output includes qcov (i.e number of genes shared between the query
            and target over number of genes on the query),
            api (average profile identity) and tapi (api x qcov)

            usage                                                                      
            -----                                                                      

            vcat utils api  [OPTIONS] -i contigs.fasta -g configs.gff -d [DBDIR]

        """,
)
@click.option(
    "-i",
    "--input",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="BLAST tabular (m8) output file",
    required=True,
)
@click.option(
    "-g",
    "--gff",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="gff file with gene cordinates of the query sequences",
    required=True,
)
@click.option(
    "-o",
    "--output",
    type=str,
    default=None,
    help="output file path",
    required=False,
)
@click.option(
    "--header",
    callback=parse_csv,
    help="columnames of the m8 file ",
    default=HEADER,
    required=False,
)
@click.option(
    "--tapif",
    type=float,
    default=0.3,
    help="assign sequences above this tapi threshold to families",
    required=False,
)
@click.option(
    "--tapio",
    type=float,
    default=0.15,
    help="assign sequences above this tapi threshold to orders",
    required=False,
)
@click.option(
    "--tapic",
    type=float,
    default=0.15,
    help="assign sequences above this tapi threshold to classes",
    required=False,
)
@click.option(
    "--tapip",
    type=float,
    default=0.15,
    help="assign sequences above this tapi threshold to phyla",
    required=False,
)
@click.option(
    "--tapik",
    type=float,
    default=0.15,
    help="assign sequences above this tapi threshold to kingdoms",
    required=False,
)
@click.option(
    "-d",
    "--dbdir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="path to the ref databases",
    required=True,
)
@click.option(
    "--topk",
    type=int,
    default=5,
    help="number of hits to select per query",
    required=False,
)
@click.option(
    "--batch",
    type=int,
    default=5000,
    help="number of records to process at a time",
    required=False,
)
@click.option(
    "--level",
    type=str,
    default=None,
    help="leave out taxonomic level (species, genus, family, class)",
    required=False,
)
@click.option(
    "--all",
    is_flag=True,
    default=False,
    help="get api for all top-k hits per query sequence. by default only outputs the besthit",
    required=False,
)
def api(input, output, header, tapif, tapio, tapic, tapip, tapik, batch, dbdir, gff, topk, level, all):
    """
    calculates the average profile identity and coverage of a query sequence
    to the genomes in the target database
    """
    CHUNK_SIZE = batch
    THRESHOLDS = {
        "family": tapif,
        "order": tapio,
        "class": tapic,
        "phylum": tapip,
        "kingdom": tapik,
    }
    file_name = os.path.basename(input)

    index = index_m8(input, kind="axi")
    tmp_files = []
    with logging_redirect_tqdm():
        for i in tqdm(range(0, len(index), CHUNK_SIZE), ncols=70, ascii=" ="):
            finput = load_chunk(
                input, index=index, recstart=i, recend=min(i + CHUNK_SIZE, len(index))
            )
            outfile = os.path.join(
                os.path.dirname(input),
                f"{os.path.splitext(file_name)[0]}_api_{min(i + CHUNK_SIZE, len(index))}.tsv",
            )
            status = axi_summary(
                finput, gff, dbdir, header, THRESHOLDS, top_k=topk, kind="api", level=level, all=all
            )
            if isinstance(status, pl.DataFrame):
                status.write_csv(outfile, separator="\t")
                logger.info(f"{outfile} updated")
                tmp_files.append(outfile)

            else:
                # exit code 1
                logger.error("error occured!")
                logger.exception(status)
                exit(1)

    logger.info("merging temporary files")
    tmp = [pl.read_csv(f, separator="\t") for f in tmp_files]
    tmp = [i for i in tmp if not i.is_empty()]
    if not output:
        outfile = os.path.join(
                os.path.dirname(input), f"{os.path.splitext(file_name)[0]}_api.tsv"
            )
    else:
        outfile = output

    if tmp:
        df = pl.concat(tmp)
        df.write_csv(outfile, separator="\t")
        logger.info(f"{outfile} updated")
    else:
        logger.info("all tables are empty")

    # Remove temporary TSV files
    for file in tmp_files:
        os.remove(file)
        Path(outfile).touch()


@utils.command(
    context_settings=dict(ignore_unknown_options=True, show_default=True),
    help="""
            subsamples nucleotide fragments from a multifasta file from a 
            given size range

            usage                                                                      
            -----                                                                      

            vcat utils fragment  [OPTIONS] -i contigs.fasta -o fragments_contig.fasta

        """,
)
@click.option(
    "-i",
    "--input",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="multi fasta file",
    required=True,
)
@click.option(
    "-o",
    "--output",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="output fasta file with sequence fragments",
    required=True,
)
@click.option(
    "--minlen",
    type=int,
    default=5000,
    help="minium fragement length",
    required=False,
)
@click.option(
    "--maxlen",
    type=int,
    default=5000,
    help="maxmium fragment length",
    required=False,
)
@click.option(
    "--overlap",
    type=int,
    default=0,
    help="number of overlapping bases between two consecutive windows",
    required=False,
)
@click.option(
    "--coverage",
    type=float,
    default=None,
    help="sample sub-sequences until this coverage threshold is met",
    required=False,
)
@click.option(
    "--circular",
    is_flag=True,
    default=False,
    help="whether the genomes are circular or not",
    required=False,
)
@click.option(
    "--man_n_prop",
    type=float,
    default=0.3,
    help="maximum proportion of Ns allowed in a fragment",
    required=False,
)
@click.option(
    "--seed",
    type=int,
    default=42,
    help="seed for the random number generator",
    required=False,
)
def fragment(**kwargs):
    """
    generates nucleotide framents from input multi fasta file (comming soon)
    """
    from vcat.frgment import split_core
    split_core(**kwargs)

@utils.command(
    context_settings=dict(ignore_unknown_options=True,  show_default=True),
    help="""
    benchmark the performance by leaving out taxa. Suppose target X belongs to taxon Y 
    we remove all query hits to taxon Y and calculate the axi to the taxa belonging to remianing hits.

    vcat contigs -i <genomes_used_to_build_db> -d <db> -o <results_dir>

    vcat utils benchmark --dbdir <path> --results <results_dir>
    """,
)
@click.option(
    "-d",
    "--dbdir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="database directory",
    required=True,
)
@click.option(
    "-r",
    "--results",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="this directory should contain vcat results",
    required=True,
)
@click.option(
    "-l",
    "--level",
    type=str,
    help="taxonomic level to leave out (species, genus, family, class)",
    required=True,
)
@click.option(
    "--tapif",
    type=float,
    default=0.30,
    help="assign sequences above this tapi threshold to families",
    required=False,
)
@click.option(
    "--tapio",
    type=float,
    default=0.15,
    help="assign sequences above this tapi threshold to orders",
    required=False,
)
@click.option(
    "--tapic",
    type=float,
    default=0.15,
    help="assign sequences above this tapi threshold to classes",
    required=False,
)
@click.option(
    "--tapip",
    type=float,
    default=0.15,
    help="assign sequences above this tapi threshold to phyla",
    required=False,
)
@click.option(
    "--tapik",
    type=float,
    default=0.15,
    help="assign sequences above this tapi threshold to kingdoms",
    required=False,
)
@click.option(
    "--taaig",
    type=float,
    default=0.49,
    help="assign sequences above this taai threshold to genera",
    required=False,
)
@click.option(
    "--taaif",
    type=float,
    default=0.3,
    help="assign sequences above this taai threshold to families",
    required=False,
)
@click.option(
    "--taaio",
    type=float,
    default=0.3,
    help="assign sequences above this taai threshold to orders",
    required=False,
)
@click.option(
    "--taaic",
    type=float,
    default=0.3,
    help="assign sequences above this taai threshold to classes",
    required=False,
)
@click.option(
    "--taaip",
    type=float,
    default=0.3,
    help="assign sequences above this taai threshold to phyla",
    required=False,
)
@click.option(
    "--taaik",
    type=float,
    default=0.3,
    help="assign sequences above this taai threshold to kingdoms",
    required=False,
)
@click.option(
    "--tanig",
    type=float,
    default=0.49,
    help="assign sequences above this tani threshold to genera",
    required=False,
)
@click.option(
    "--tanis",
    type=float,
    default=0.81,
    help="assign sequences above this tani threshold to species",
    required=False,
)
@click.option(
    "--batch",
    type=int,
    default=5000,
    help="number of records to process at a time",
    required=False,
)
@click.option(
    "-n",
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="Test execution.",
)

@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def benchmark(dbdir, results, batch, level, snakemake_args, **kwargs):
    """
    benchmark the performance by leaving out taxa. Suppose target X belongs to taxon Y 
    we remove all query hits to taxon Y and calculate the axi to the taxa belonging to remianing hits.

    vcat contigs -i <genomes_used_to_build_db> -d <db> -o <results_dir>

    vcat utils benchmark --dbdir <path> --results <results_dir>
    """
    pass
    logger.info(f"vcat version: {__version__}")
    taai_params = ""
    tapi_params = ""
    tani_params = ""
    for k, v in kwargs.items():
        if k.startswith("taai"):
            taai_params += f" --{k} {v}"
        elif k.startswith("tapi"):
            tapi_params += f" --{k} {v}"
        elif k.startswith("tani"):
            tani_params += f" --{k} {v}"
        

    jobs = kwargs.get("jobs", 4)
    cmd_parts = [
        "snakemake",
        f"--snakefile {get_snakefile('./pipeline/rules/benchmark.smk')}",
        f"--jobs {jobs}",
        "--rerun-incomplete",
        "--scheduler greedy",
        "--show-failed-logs",
        "--groups group1=1",
        "--config",
        f"database_dir='{dbdir}'",
        f"results='{results}'",
        f"batch='{batch}'",
        f"ani='{tani_params}'",
        f"aai='{taai_params}'",
        f"api='{tapi_params}'",
        f"level='{level}'",
    ]

    if snakemake_args:
        if snakemake_args[0].startswith("-"):
            cmd_parts.extend(snakemake_args)
        else:
            cmd_parts.append("--")
            cmd_parts.extend(snakemake_args)

    cmd = " ".join(cmd_parts)

    logger.info(cmd)
    logger.debug("Executing: %s" % cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logger.critical(e)
        exit(1)

cli.add_command(utils)

if __name__ == "__main__":
    cli()
