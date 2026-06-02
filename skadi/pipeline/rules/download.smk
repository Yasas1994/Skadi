from pathlib import Path
import os

# Set the input file and output directory as global variables
configfile: "config.yaml"
# these values are incuded in the snakefile
DBDIR = config["database_dir"]
DBURL = config["dburl"]
DBVER = config["dbversion"]

rule all:
    input:
        f"{DBDIR}/{DBVER}/.done"

# Download the pre-build database
rule download_db:
    output:
        f"{DBDIR}/{DBVER}.tar.gz"
    params:
        url =DBURL,
	    outdir = DBDIR,
	    outfile = f"{DBVER}.tar.gz"
    shell:
        """
        wget {params.url} -O {output}
        """

# Download the pre-build database
rule decompress:
    input:
        f"{DBDIR}/{DBVER}.tar.gz"
    output:
        f"{DBDIR}/{DBVER}/.done"
    params:
        outdir = DBDIR
    shell:
        """
        tar -xf {input} -C {params.outdir} && touch {output}
        rm {input}
        """
