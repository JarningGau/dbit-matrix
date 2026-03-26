# Developer Architecture

This page is the developer-facing overview of the repository.

## Workflow Families

The repo provides three orchestration entrypoints:

- `scripts/make_cmd.py` for `TAPS`
- `scripts-emseq/make_cmd.py` for `EMSeq`
- `scripts-rna/make_cmd.py` for `RNA`

Stage implementations stay in thin single-step scripts:

- `scripts/*.py` for TAPS and shared utilities
- `scripts-emseq/*.py` for EMSeq-specific stages
- `scripts-rna/*.py` for RNA-specific stages

## Main Pipelines

TAPS and EMSeq mainline:

`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> saturation -> summary`

RNA current scope:

`demux_extract_bc -> align`

For EMSeq, `--stage all` covers the fixed mainline above and does not include optional experimental stages.

## Execution Model

- all workflow parameters should live in `workflow/*.json`
- all stages support `--dry-run`
- `local` generates shell scripts under `work/<sample>/commands/`
- `slurm` generates `.sbatch` files under `work/<sample>/commands/`
- Slurm resource settings are stage-specific rather than shared globally

## Shared Output Layout

Common output root:

- `work/<sample>/`

Frequently used subdirectories:

- `commands/`
- `demux/`
- `align_shards/`
- `pooled/`
- `split_bams/`
- `coverage/`
- `qc/`
- `summary/`
- `solo/` for RNA

## Workflow-Specific Notes

- EMSeq `mbias` is implemented by `scripts-emseq/mbias.py` and does not use the same chemistry logic as TAPS `scripts/mbias.py`
- RNA currently produces STARsolo matrix outputs rather than methylation coverage outputs
