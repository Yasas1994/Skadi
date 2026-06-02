## Installation


currently you can install the development version of skadi with following steps.


#### Installing with conda/mamba
clone the git repositoty and use the environment.yaml file to create a conda environemnt.
alternatively you can also use pixi. Then, install the skadi python package into the conda 
environment you just created.

```bash
git clone https://github.com/Yasas1994/skadi.git
cd skadi
mamba create -f environment.yml
mamba activate skadi
pip install .
skadi --help
```



#### Running skadi with singularity
clone the git reposity and used the Apptainer definition file to build a singularity container.

```bash
git clone https://github.com/Yasas1994/skadi.git
cd skadi
apptainer build skadi.sif Apptainer
apptainer run skadi.sif skadi --help
```


if everything goes soothly, you should see skadi help on the the terminal.

```text
Usage: skadi [OPTIONS] COMMAND [ARGS]...

  skadi: a command-line tool-kit for adding ICTV taxonomy annotations to virus
  contigs, mapping reads to virus genomes and much more.
  (https://github.com/Yasas1994/skadi)

Options:
  --version   Show the version and exit.
  -h, --help  Show this message and exit.

Commands:
  contigs    run contig annotation workflow
  downloaddb download pre-built reference databases
  preparedb  download sequences and build reference databases
  reads      run read annotation workflow
  utils      tool chain for calculating ani, aai and visualizations
```


#### Downloading pre-built databases

You can download pre-built databases instead of building from scratch. 
Currently, msl39v and masl40v1 are available to download

```bash

skadi downloaddb --dbversion masl40v1 -d <path-to-save-the-database> --cores 1
```

optionally, you can also build it yourself with the following command.

```bash
skadi preparedb -d <path-to-save-the-database>
```
