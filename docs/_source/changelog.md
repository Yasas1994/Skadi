# Changelog

All notable changes to this project will be documented in this file.

## 0.0.4 — 2026-06-05

### Added
- Per-family threshold tuning for ANI, AAI, and API scores.
- `benchmark` subcommand for threshold evaluation.
- Cluster deployment guide and SLURM job templates.
- Apptainer/Singularity container support with editable installs.
- `cluster-setup` documentation for HPC deployment.

### Fixed
- `--max-members None` YAML null parsing issue.
- Snakemake `{config.get(...)}` formatting in shell commands.
- `load_chunk` index key mismatch in post-processing.
- `rank_taxid` Int32/Int64 schema error.
- `taai`/`tapi` column rename consistency.
- Container read-only script path handling.

### Changed
- Updated database versions: msl39v4, msl40v2, msl41v1.

## 0.0.3 — 2025-03-15

### Added
- msl40v2 is now available to download.
- Added read annotation workflow (`skadi reads`).

## 0.0.2 — 2025-01-20

### Added
- `downloaddb` command to download pre-built databases (`msl39v4` and `msl40v1`).
- Apptainer definition file for container builds.
- Subroutines to clean up temporary files generated during database build.
