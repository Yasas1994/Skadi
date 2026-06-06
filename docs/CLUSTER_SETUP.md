# SKADI Cluster Deployment Guide

This guide describes how to deploy and run SKADI on HPC clusters using Apptainer/Singularity containers and workload managers like Slurm.

## Overview

Many HPC clusters restrict software installation on login nodes and use workload managers (Slurm, PBS, LSF) for job scheduling. This guide provides a **container-based deployment** that works on any cluster with:

- Apptainer or Singularity installed
- A workload manager (Slurm examples provided)
- A shared filesystem (NFS, Lustre, BeeGFS, etc.)

## Architecture

```
┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
│   Local Machine │────▶│  Cluster Login  │────▶│  Compute Nodes  │
│  (dev + rsync)  │     │  (ssh + submit) │     │  (Container)    │
└─────────────────┘     └─────────────────┘     └─────────────────┘
                                                        │
                              ┌─────────────────────────┘
                              ▼
                    ┌─────────────────┐
                    │  skadi.sif      │  (Apptainer container)
                    │  + bind mounts  │
                    │  /skadi (code)  │
                    │  /data/db       │
                    │  /data/input    │
                    │  /data/output   │
                    └─────────────────┘
```

## Prerequisites

### On the cluster

- **Apptainer** (≥1.0) or **Singularity** (≥3.0)
- **Slurm** (or other workload manager)
- **Shared filesystem** accessible from all nodes
- **Internet access** on compute nodes (for initial database build) OR pre-downloaded data

### On your local machine

- Apptainer/Singularity (for local container testing, optional)
- `rsync` or `scp` for file transfer
- SSH access to the cluster

## Quick Start

```bash
# 1. SSH to cluster login node
ssh your-cluster

# 2. Navigate to SKADI directory
cd /path/to/skadi

# 3. Build database (one-time, ~30 min, 128GB RAM recommended)
sbatch slurm/preparedb.slurm

# 4. Run contig annotation
sbatch slurm/contigs.slurm

# 5. Check job status
squeue -u $USER
```

## Step-by-Step Setup

### 1. Choose a Base Directory

Select a path on the shared filesystem:

```bash
export SKADI_DIR=/path/to/your/skadi
mkdir -p $SKADI_DIR/{code,db,data,output,logs,slurm}
```

### 2. Create the Container Definition

Create `singularity/skadi.def`:

```singularity
Bootstrap: docker
From: condaforge/miniforge3:latest

%files
    ./environment.yml /app/environment.yml
    ./setup.py /app/setup.py
    ./pyproject.toml /app/pyproject.toml
    ./README.md /app/README.md
    ./MANIFEST.in /app/MANIFEST.in
    ./LICENSE /app/LICENSE
    ./scripts /app/scripts
    ./skadi /app/skadi

%post
    apt-get update && apt-get install -y --no-install-recommends \
        git wget curl ca-certificates && rm -rf /var/lib/apt/lists/*
    conda env create -n skadi --file /app/environment.yml
    conda clean --all --yes
    /opt/conda/envs/skadi/bin/pip install --root-user-action=ignore /app

%environment
    export PATH="/opt/conda/envs/skadi/bin:$PATH"
    export CONDA_DEFAULT_ENV=skadi

%labels
    Author "Your Name"
    Version "0.0.4"
    Description "SKADI pipeline container"
```

### 3. Build the Container

**Important**: Build on a compute node, not the login node, to avoid memory issues with `mksquashfs`.

Create `slurm/build_container.slurm`:

```bash
#!/bin/bash
#SBATCH --job-name=skadi_build
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=logs/build_%j.log

set -e
cd $SKADI_DIR/code
rm -f $SKADI_DIR/skadi.sif
apptainer build $SKADI_DIR/skadi.sif singularity/skadi.def
echo "Build complete: $(ls -lh $SKADI_DIR/skadi.sif)"
```

Submit:
```bash
sbatch slurm/build_container.slurm
```

### 4. Create Cluster-Specific Configuration

Create `skadi/pipeline/config_cluster.yaml`:

```yaml
database_dir: /path/to/your/skadi/db
# Keep all other settings from config.yaml
downloads:
  - latest: true
    link: https://ckan.fdm.uni-greifswald.de/.../msl41v1.tar.gz
    name: msl41v1
profile_cluster:
  min_seq_id: 0.0
  coverage: 0.5
  cov_mode: 0
  sensitivity: 7.5
  cluster_mode: 0
  threads: 8
lca:
  fraction: 0.6
profile_msa:
  qid: 0.0
  qsc: -20.0
  cov: 0.0
  max_seq_id: 0.9
  filter_msa: 1
  e_profile: 0.001
cluster_filter:
  min_members: 1
  max_members: null
```

### 5. Create SLURM Scripts

#### `preparedb.slurm`

```bash
#!/bin/bash
#SBATCH --job-name=skadi_preparedb
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G        # REQUIRED: diamond makedb needs >64GB
#SBATCH --time=02:00:00
#SBATCH --output=logs/preparedb_%j.log

set -e
echo "Running on: $(hostname)"
echo "Started at: $(date)"

SIF=$SKADI_DIR/skadi.sif
CODE=$SKADI_DIR/code
DB=$SKADI_DIR/db

# Job-specific editable install
INSTALL_DIR=$SKADI_DIR/.local_${SLURM_JOB_ID}
mkdir -p $INSTALL_DIR/lib/python3.11/site-packages
export PYTHONPATH="$INSTALL_DIR/lib/python3.11/site-packages:$PYTHONPATH"

mkdir -p $DB $SKADI_DIR/logs

apptainer exec \
  --bind $CODE:/skadi \
  --bind $DB:/data/db \
  --bind $INSTALL_DIR:/home/$USER/.local \
  $SIF bash -c "
    cd /skadi
    pip install --target=/home/$USER/.local/lib/python3.11/site-packages --no-deps -e . >/dev/null 2>&1
    skadi preparedb -d /data/db --jobs 16 \
      --configfile /skadi/skadi/pipeline/config_cluster.yaml
  "

echo "Finished at: $(date)"
```

#### `contigs.slurm`

```bash
#!/bin/bash
#SBATCH --job-name=skadi_contigs
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G         # REQUIRED: mmseqs profile search needs >32GB
#SBATCH --time=04:00:00
#SBATCH --output=logs/contigs_%j.log

set -e
echo "Running on: $(hostname)"
echo "Started at: $(date)"

SIF=$SKADI_DIR/skadi.sif
CODE=$SKADI_DIR/code
DB=$SKADI_DIR/db
INPUT=$SKADI_DIR/data/your_contigs.fna
OUTPUT=$SKADI_DIR/output/your_results

INSTALL_DIR=$SKADI_DIR/.local_${SLURM_JOB_ID}
mkdir -p $INSTALL_DIR/lib/python3.11/site-packages
export PYTHONPATH="$INSTALL_DIR/lib/python3.11/site-packages:$PYTHONPATH"

mkdir -p $SKADI_DIR/logs $SKADI_DIR/output

apptainer exec \
  --bind $CODE:/skadi \
  --bind $DB:/data/db \
  --bind $(dirname $INPUT):/data/input \
  --bind $(dirname $OUTPUT):/data/output \
  --bind $INSTALL_DIR:/home/$USER/.local \
  $SIF bash -c "
    cd /skadi
    pip install --target=/home/$USER/.local/lib/python3.11/site-packages --no-deps -e . >/dev/null 2>&1
    skadi contigs -i /data/input/$(basename $INPUT) -o /data/output/$(basename $OUTPUT) \
      -d /data/db --jobs 8 \
      --configfile /skadi/skadi/pipeline/config_cluster.yaml
  "

echo "Finished at: $(date)"
```

#### `reads.slurm`

```bash
#!/bin/bash
#SBATCH --job-name=skadi_reads
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/reads_%j.log

set -e
SIF=$SKADI_DIR/skadi.sif
CODE=$SKADI_DIR/code
DB=$SKADI_DIR/db
INPUT=$SKADI_DIR/data/reads_1.fastq
INPUT2=${2:-}  # Optional paired-end
OUTPUT=$SKADI_DIR/output/reads_results

INSTALL_DIR=$SKADI_DIR/.local_${SLURM_JOB_ID}
mkdir -p $INSTALL_DIR/lib/python3.11/site-packages
export PYTHONPATH="$INSTALL_DIR/lib/python3.11/site-packages:$PYTHONPATH"

# Build input bind and args
INPUT_BIND="--bind $(dirname $INPUT):/data/input"
INPUT_ARG="-i /data/input/$(basename $INPUT)"
if [ -n "$INPUT2" ]; then
    INPUT_BIND="$INPUT_BIND --bind $(dirname $INPUT2):/data/input2"
    INPUT_ARG="$INPUT_ARG -i2 /data/input2/$(basename $INPUT2)"
fi

apptainer exec \
  --bind $CODE:/skadi \
  --bind $DB:/data/db \
  $INPUT_BIND \
  --bind $(dirname $OUTPUT):/data/output \
  --bind $INSTALL_DIR:/home/$USER/.local \
  $SIF bash -c "
    cd /skadi
    pip install --target=/home/$USER/.local/lib/python3.11/site-packages --no-deps -e . >/dev/null 2>&1
    skadi reads $INPUT_ARG -o /data/output/$(basename $OUTPUT) \
      -d /data/db --jobs 8 \
      --configfile /skadi/skadi/pipeline/config_cluster.yaml
  "
```

### 6. Sync Code to Cluster

From your local machine:

```bash
rsync -av --exclude='.git' --exclude='.snakemake' --exclude='__pycache__' \
  --exclude='*.egg-info' ./ your-cluster:$SKADI_DIR/code/
```

### 7. Build Database

```bash
ssh your-cluster
cd $SKADI_DIR
sbatch slurm/preparedb.slurm
```

Monitor with `squeue -u $USER` and check `logs/preparedb_*.log`.

### 8. Run Analysis

```bash
# Contigs
sbatch slurm/contigs.slurm

# Reads (single-end)
sbatch slurm/reads.slurm /path/to/reads.fastq /path/to/output

# Reads (paired-end)
sbatch slurm/reads.slurm /path/to/R1.fastq /path/to/output /path/to/R2.fastq
```

## Resource Requirements

| Step | Min Memory | Recommended | Notes |
|------|-----------|-------------|-------|
| `diamond makedb` | 64GB | **128GB** | Taxonomy loading is memory-intensive |
| `mmseqs profile search` | 32GB | **64GB** | Prefilter needs substantial RAM |
| `mmseqs cluster` | 16GB | 32GB | For profile clustering |
| `skadi contigs` (total) | 32GB | **64GB** | Covers all pipeline steps |
| `skadi reads` (total) | 16GB | 32GB | BBMap + alignment steps |
| `skadi preparedb` | 64GB | **128GB** | Full database build |

## Updating Code

### Python package changes (`skadi/`)

The job-specific pip install in SLURM scripts picks up changes automatically. Just sync and re-submit:

```bash
rsync -av ./skadi/ your-cluster:$SKADI_DIR/code/skadi/
```

### Script changes (`scripts/`)

Scripts installed to `/opt/conda/envs/skadi/bin/` inside the container are **baked in at build time**. To update:

**Option A**: Rebuild container
```bash
sbatch slurm/build_container.slurm
```

**Option B**: Call from bind mount (modify Snakefile)
```python
# In your Snakefile, change:
#   postprocess.py {params} ...
# To:
#   python3 /skadi/scripts/postprocess.py {params} ...
```

### Configuration changes

Edit `skadi/pipeline/config_cluster.yaml` and sync. No rebuild needed.

## Troubleshooting

### OOM Kill (exit 137)

Increase `--mem` in your SLURM script. Common culprits:
- `diamond makedb` with taxonomy: needs 128GB
- `mmseqs search` prefilter: needs 64GB

### "Prefilter died" / MMseqs errors

MMseqs ran out of memory. Either:
- Increase SLURM `--mem`
- Reduce `--split-memory-limit` in the Snakefile
- Use fewer threads (`--jobs`)

### Container build fails on login node

Always build via `sbatch` on a compute node. The `mksquashfs` step needs significant memory.

### "No such file" inside container

Check your bind mounts. Paths inside the container (`/skadi`, `/data/db`) must match what the Snakefile expects.

### Slow filesystem performance

If your cluster uses NFS or Lustre, consider:
- Setting `SNAKEMAKE_PROFILE` to use local temp directories
- Adding `--shadow-prefix /tmp` to Snakemake calls
- Increasing `--latency-wait` for slow filesystems

## Adapting to Other Workload Managers

### PBS/Torque

Replace `sbatch` with `qsub` and `#SBATCH` with `#PBS`:

```bash
#PBS -N skadi_preparedb
#PBS -l nodes=1:ppn=16
#PBS -l mem=128gb
#PBS -l walltime=02:00:00
#PBS -o logs/preparedb.log
```

### LSF

```bash
#BSUB -J skadi_preparedb
#BSUB -n 16
#BSUB -M 128000
#BSUB -W 02:00
#BSUB -o logs/preparedb.log
```

## Advanced Topics

### GPU Support

If your cluster has GPU nodes and you need GPU-accelerated tools:

```bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

apptainer exec --nv \  # Enable NVIDIA support
  --bind ... \
  $SIF bash -c "..."
```

### Multi-Node Execution

SKADI's Snakemake pipeline can use cluster profiles for multi-node execution:

```bash
skadi preparedb -d /data/db --jobs 64 \
  --profile /path/to/slurm-profile
```

See [Snakemake cluster profiles](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) for details.

### Caching Container Images

On clusters with local node storage, copy the SIF to `$TMPDIR` or `/tmp` at job start:

```bash
# In your SLURM script
LOCAL_SIF=/tmp/skadi_$$.sif
cp $SKADI_DIR/skadi.sif $LOCAL_SIF
# Use $LOCAL_SIF in apptainer exec
```

## References

- [Apptainer Documentation](https://apptainer.org/docs/)
- [Snakemake Cluster Execution](https://snakemake.readthedocs.io/en/stable/executing/cluster.html)
- [Slurm Documentation](https://slurm.schedmd.com/documentation.html)
