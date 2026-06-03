import click
import subprocess
import multiprocessing
import os
import re
import sys
import yaml
import polars as pl
from pathlib import Path
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
from skadi.utils import ani_summary, axi_summary, index_m8, load_chunk
from .color_logger import logger
from importlib.metadata import version, PackageNotFoundError


# Characters that are dangerous in shell arguments or paths
_UNSAFE_PATH_CHARS = re.compile(r"[;&|`$(){}\[\]\\\n\r\x00]")


def _validate_path(path: str | os.PathLike, name: str = "path") -> str:
    """Validate a user-provided path for suspicious characters.

    Args:
        path: The path to validate.
        name: Human-readable name for error messages.

    Returns:
        The path as a string.

    Raises:
        click.BadParameter: If the path contains unsafe characters.
    """
    path_str = str(path)
    if _UNSAFE_PATH_CHARS.search(path_str):
        raise click.BadParameter(
            f"{name} contains unsafe characters: {path_str!r}"
        )
    return path_str


# Define the directory containing the pipeline files
PIPELINE_DIR = os.path.join(os.path.dirname(__file__), "./pipeline")
CONFIG = os.path.join(PIPELINE_DIR, "config.yaml")
HEADER = "query,target,theader,fident,qlen,tlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,taxid,taxname,taxlineage"
CONFIG_CONTENT = yaml.safe_load(open(CONFIG, "r"))

__version__ = ""
try:
    __version__ = version("skadi")
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


def get_snakefile(file=f"{PIPELINE_DIR}/Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf


def update_config(config_path, data: dict):
    import yaml

    with open(config_path, "r") as file:
        config = yaml.safe_load(file)

    for k, v in data.items():
        config[k] = v

    with open(config_path, "w") as file:
        yaml.safe_dump(config, file, default_flow_style=False)


def _run_command(cmd: list[str]) -> None:
    """Run a command safely without shell=True.

    Args:
        cmd: List of command arguments (e.g., ["snakemake", "--jobs", "4"]).

    Raises:
        SystemExit: If the command returns a non-zero exit code.
    """
    logger.debug("Executing: %s", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.critical("Command failed: %s", e)
        sys.exit(1)
    except FileNotFoundError as e:
        logger.critical("Command not found: %s", e)
        sys.exit(1)


def _build_snakemake_cmd(
    snakefile: str,
    jobs: int,
    configfile: str,
    config_overrides: dict[str, str],
    extra_args: tuple[str, ...] = (),
    dryrun: bool = False,
) -> list[str]:
    """Build a snakemake command as a list for safe execution.

    Args:
        snakefile: Path to the Snakefile.
        jobs: Number of parallel jobs.
        configfile: Path to the config YAML file.
        config_overrides: Dict of --config key=value pairs.
        extra_args: Additional snakemake arguments.
        dryrun: Whether to add --dry-run.

    Returns:
        Command as a list of strings.
    """
    cmd = [
        "snakemake",
        "--snakefile", snakefile,
        "--jobs", str(jobs),
        "--rerun-incomplete",
        "--configfile", configfile,
        "--scheduler", "greedy",
        "--show-failed-logs",
        "--groups", "group1=1",
        "--config",
    ]
    for key, value in config_overrides.items():
        cmd.append(f"{key}={value}")

    if dryrun:
        cmd.append("--dry-run")

    if extra_args:
        if extra_args[0].startswith("-"):
            cmd.extend(extra_args)
        else:
            cmd.append("--")
            cmd.extend(extra_args)

    return cmd


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context
def cli(obj):
    r"""
    SKADI: Sequence-based Knowledgebase for Annotation, Detection, and Identification
    (https://github.com/Yasas1994/skadi)"""
    pass


@cli.command(
    context_settings=dict(ignore_unknown_options=True, show_default=True),
    short_help="run contig annotation workflow",
    help="""
    SKADI: Sequence-based Knowledgebase for Annotation, Detection, and Identification

    """,
)
@click.option(
    "-i",
    "--input",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="input file/s to run skadi on",
    required=True,
)
@click.option(
    "-o",
    "--output",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="dir to store skadi results",
    required=True,
)
@click.option(
    "-d",
    "--database",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="dir to skadi database",
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
def contigs(input, output, database, jobs, batch, snakemake_args, dryrun, **kwargs):
    """
    Runs SKADI pipeline on contigs

    Most snakemake arguments can be appended to the command for more info see 'snakemake --help'
    """
    _validate_path(input, name="input")
    _validate_path(output, name="output")
    if database:
        _validate_path(database, name="database")

    logger.info(f"skadi version: {__version__}")
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

    nuc_search = kwargs.get("nuc_search", None)

    config_overrides = {
        "database_dir": db_dir,
        "sample": input,
        "output_dir": output,
        "api": tapi_params,
        "aai": taai_params,
        "ani": tani_params,
        "batch": str(batch),
    }
    if nuc_search is not None:
        config_overrides["nuc_search"] = nuc_search

    cmd = _build_snakemake_cmd(
        snakefile=get_snakefile("./pipeline/contig_annotation.smk"),
        jobs=jobs,
        configfile=CONFIG,
        config_overrides=config_overrides,
        extra_args=snakemake_args,
        dryrun=dryrun,
    )

    logger.info(" ".join(cmd))
    _run_command(cmd)


@cli.command(
    context_settings=dict(ignore_unknown_options=True, show_default=True),
    short_help="run read annotation workflow",
    help="""
    SKADI: Sequence-based Knowledgebase for Annotation, Detection, and Identification

    usage (paired-end): skadi reads [OPTIONS] -i1 pair1.fastq -i2 pair2.fastq -o mapping_results.tsv

    usage (single-end): skadi reads [OPTIONS] -i1 pair1.fastq -o mapping_results.tsv
    """,
)
@click.option(
    "-in",
    "--input",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="input read file1 to run skadi on",
    required=True,
)
@click.option(
    "-in2",
    "--input2",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="input read file2 (paired-end) to run skadi on",
    required=False,
)
@click.option(
    "-o",
    "--output",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="dir to store skadi results",
    required=True,
)
@click.option(
    "-d",
    "--database",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="dir to skadi database",
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
    Runs SKADI pipeline on reads

    Most snakemake arguments can be appended to the command for more info see 'snakemake --help'
    """
    _validate_path(input, name="input")
    if input2:
        _validate_path(input2, name="input2")
    _validate_path(output, name="output")
    if database:
        _validate_path(database, name="database")

    logger.info(f"skadi version: {__version__}")
    conf = load_configfile(CONFIG)
    db_dir = database or conf["database_dir"]
    logger.info(f"database version: {Path(db_dir).name}")

    sample_str = f"in={input}, in2={input2}" if input2 else f"in={input}"

    config_overrides = {
        "database_dir": db_dir,
        "sample": sample_str,
        "output_dir": output,
        "bbmap_args": bbmap_args,
        "pileup_args": pileup_args,
        "summary_args": summary_args,
    }

    cmd = _build_snakemake_cmd(
        snakefile=get_snakefile("./pipeline/read_annotation.smk"),
        jobs=jobs,
        configfile=CONFIG,
        config_overrides=config_overrides,
        extra_args=snakemake_args,
        dryrun=dryrun,
    )

    if profile:
        # Insert --profile after snakemake but before other args
        cmd.insert(2, "--profile")
        cmd.insert(3, profile)

    logger.debug("Executing: %s", " ".join(cmd))
    _run_command(cmd)


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
    """Executes a snakemake workflow to download and build the databases"""
    _validate_path(db_dir, name="db-dir")

    logger.info("Building taxdb")
    cmd = _build_snakemake_cmd(
        snakefile=get_snakefile("./pipeline/rules/taxdump.smk"),
        jobs=jobs,
        configfile=CONFIG,
        config_overrides={"database_dir": db_dir},
        extra_args=snakemake_args,
    )
    _run_command(cmd)

    logger.info("Building mmseqs databases")
    cmd = _build_snakemake_cmd(
        snakefile=get_snakefile("./pipeline/rules/createdb.smk"),
        jobs=jobs,
        configfile=CONFIG,
        config_overrides={"database_dir": db_dir},
        extra_args=snakemake_args,
    )
    _run_command(cmd)

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
    _validate_path(db_dir, name="db-dir")

    logger.info(f"downloading {dbversion} database from remote server")
    cmd = _build_snakemake_cmd(
        snakefile=get_snakefile("./pipeline/rules/download.smk"),
        jobs=1,
        configfile=CONFIG,
        config_overrides={
            "database_dir": db_dir,
            "dbversion": dbversion,
            "dburl": default_url,
        },
        extra_args=snakemake_args,
    )
    _run_command(cmd)

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

            skadi utils ani [OPTIONS] -i contigs.fasta

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
    _validate_path(input, name="input")
    if output:
        _validate_path(os.path.dirname(output) or ".", name="output directory")
    if dbdir:
        _validate_path(dbdir, name="dbdir")

    if (level or dbdir) and not (level and dbdir):
        logger.error(f"{level} and {dbdir} are mutually inclusive")
        sys.exit(1)

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
            try:
                status = ani_summary(finput, all=all, header=header, level=level, dbdir=dbdir)
            except Exception as e:
                logger.error("error occured!")
                logger.exception(e)
                sys.exit(1)

            if isinstance(status, pl.DataFrame):
                status.write_csv(outfile, separator="\t")
                logger.info(f"{outfile} updated")
                tmp_files.append(outfile)
            else:
                logger.error("ani_summary returned unexpected type: %s", type(status))
                sys.exit(1)

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

            skadi utils aai  [OPTIONS] -i contigs.fasta -g configs.gff -d [DBDIR]

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
    _validate_path(input, name="input")
    _validate_path(gff, name="gff")
    _validate_path(dbdir, name="dbdir")
    if output:
        _validate_path(os.path.dirname(output) or ".", name="output directory")

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
            try:
                status = axi_summary(
                    finput, gff, dbdir, header, THRESHOLDS, top_k=topk, kind="aai", all=all, level=level
                )
            except Exception as e:
                logger.error("error occured!")
                logger.exception(e)
                sys.exit(1)

            if isinstance(status, pl.DataFrame):
                status.write_csv(outfile, separator="\t")
                logger.info(f"{outfile} updated")
                tmp_files.append(outfile)
            else:
                logger.error("axi_summary returned unexpected type: %s", type(status))
                sys.exit(1)

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

            skadi utils api  [OPTIONS] -i contigs.fasta -g configs.gff -d [DBDIR]

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
    _validate_path(input, name="input")
    _validate_path(gff, name="gff")
    _validate_path(dbdir, name="dbdir")
    if output:
        _validate_path(os.path.dirname(output) or ".", name="output directory")

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
            try:
                status = axi_summary(
                    finput, gff, dbdir, header, THRESHOLDS, top_k=topk, kind="api", level=level, all=all
                )
            except Exception as e:
                logger.error("error occured!")
                logger.exception(e)
                sys.exit(1)

            if isinstance(status, pl.DataFrame):
                status.write_csv(outfile, separator="\t")
                logger.info(f"{outfile} updated")
                tmp_files.append(outfile)
            else:
                logger.error("axi_summary returned unexpected type: %s", type(status))
                sys.exit(1)

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

            skadi utils fragment  [OPTIONS] -i contigs.fasta -o fragments_contig.fasta

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
    "--max_n_prop",
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
    _validate_path(kwargs.get("input"), name="input")
    _validate_path(kwargs.get("output"), name="output")
    from skadi.frgment import split_core
    split_core(**kwargs)


@utils.command(
    context_settings=dict(ignore_unknown_options=True, show_default=True),
    help="""
    benchmark the performance by leaving out taxa. Suppose target X belongs to taxon Y 
    we remove all query hits to taxon Y and calculate the axi to the taxa belonging to remianing hits.

    skadi contigs -i <genomes_used_to_build_db> -d <db> -o <results_dir>

    skadi utils benchmark --dbdir <path> --results <results_dir>
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
    help="this directory should contain skadi results",
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
    "--thresholds-file",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    default=None,
    help="JSON file with per-group thresholds from compute_pairwise_matrix.py",
)
@click.option(
    "--batch",
    type=int,
    default=5000,
    help="number of records to process at a time",
    required=False,
)
@click.option(
    "--method",
    type=click.Choice(["cascade", "consensus"], case_sensitive=False),
    default="cascade",
    help="Taxonomy assignment method for postprocessing.",
)
@click.option(
    "--min-confidence",
    type=float,
    default=0.5,
    help="Minimum confidence for consensus assignment.",
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
def benchmark(dbdir, results, batch, level, method, min_confidence, snakemake_args, **kwargs):
    """
    benchmark the performance by leaving out taxa. Suppose target X belongs to taxon Y 
    we remove all query hits to taxon Y and calculate the axi to the taxa belonging to remianing hits.

    skadi contigs -i <genomes_used_to_build_db> -d <db> -o <results_dir>

    skadi utils benchmark --dbdir <path> --results <results_dir>
    """
    _validate_path(dbdir, name="dbdir")
    _validate_path(results, name="results")

    logger.info(f"skadi version: {__version__}")
    taai_params = ""
    tapi_params = ""
    tani_params = ""
    thresholds_file_param = ""
    for k, v in kwargs.items():
        if k.startswith("taai"):
            taai_params += f" --{k} {v}"
        elif k.startswith("tapi"):
            tapi_params += f" --{k} {v}"
        elif k.startswith("tani"):
            tani_params += f" --{k} {v}"
        elif k == "thresholds_file" and v:
            thresholds_file_param = f" --thresholds-file {v}"

    jobs = kwargs.get("jobs", 4)
    config_overrides = {
        "database_dir": dbdir,
        "results": results,
        "batch": str(batch),
        "ani": tani_params + thresholds_file_param,
        "aai": taai_params,
        "api": tapi_params,
        "level": level,
        "method": method,
        "min_confidence": str(min_confidence),
    }

    cmd = _build_snakemake_cmd(
        snakefile=get_snakefile("./pipeline/rules/benchmark.smk"),
        jobs=jobs,
        configfile=CONFIG,
        config_overrides=config_overrides,
        extra_args=snakemake_args,
        dryrun=kwargs.get("dryrun", False),
    )

    logger.info(" ".join(cmd))
    _run_command(cmd)


cli.add_command(utils)

if __name__ == "__main__":
    cli()
