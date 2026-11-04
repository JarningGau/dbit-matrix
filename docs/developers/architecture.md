# Developer Architecture

This page is the developer-facing overview of the repository.

## Workflow Families

The repo provides four orchestration entrypoints:

- `scripts/make_cmd.py` for `TAPS`
- `scripts-v2/make_cmd.py` for `TAPS v2` (forked demux and summary; reuses most `scripts/*.py`)
- `scripts-emseq/make_cmd.py` for `EMSeq`
- `scripts-rna/make_cmd.py` for `RNA`

Stage implementations stay in thin single-step scripts:

- `scripts/*.py` for TAPS and shared utilities
- `scripts-v2/*.py` for TAPS v2 demux and summary
- `scripts-emseq/*.py` for EMSeq-specific stages
- `scripts-rna/*.py` for RNA-specific stages

## Main Pipelines

TAPS, TAPS v2, and EMSeq share the same mainline stage order. Full order and optional stages: [stage contracts](contracts.md).

TAPS / TAPS v2 may insert `hmc_spike` after `saturation` and before `summary` when `hmc_spike.enabled` is true. EMSeq does not implement `hmc_spike`.

RNA current scope:

`demux_extract_bc -> align`

For EMSeq, `--stage all` covers the fixed mainline and does not include optional experimental stages.

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
- `coverage/`
- `qc/`
- `summary/`
- `solo/` for RNA

## Workflow-Specific Notes

- EMSeq `mbias` is implemented by `scripts-emseq/mbias.py` and does not use the same chemistry logic as TAPS `scripts/mbias.py`
- EMSeq `call` is pooled pysam (`scripts-emseq/methy_caller.py` / `scripts-emseq/call.py`); biscuit remains an align-stage tool
- RNA currently produces STARsolo matrix outputs rather than methylation coverage outputs
