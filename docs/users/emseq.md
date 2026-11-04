# EMSeq User Guide

Use this guide for the EMSeq workflow.

## What It Runs

The mainline stages are `fastp_split` through `summary`. Full order and optional stages: [stage contracts](../developers/contracts.md).

`--stage all` runs the mainline only (no `hmc_spike`).

## Before You Start

Prepare:

- `R1` and `R2` FASTQ files
- `barcode1_whitelist`
- `barcode2_whitelist`
- `split_barcodes` (CB whitelist for `call`)
- `biscuit_reference`
- `call_reference_file`
- `call_jobs`
- optional spike-in reference files

Start from:

- config template: [`workflow/dbit_emseq_test.json`](../../workflow/dbit_emseq_test.json)
- entry script: `scripts-emseq/make_cmd.py`

## Minimal Setup

Install the environment:

```bash
pixi install
```

Copy and edit the example workflow JSON. At minimum, set:

- `sample_id`
- `r1`
- `r2`
- `barcode1_whitelist`
- `barcode2_whitelist`
- `linker_bc`
- `insert_left`
- `split_barcodes` (CB whitelist for `call`)
- `biscuit_reference`
- `call_reference_file`
- `call_jobs`

## Dry-Run

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage all \
  --runner local \
  --dry-run
```

If you use Slurm, change `--runner local` to `--runner slurm`.

## Run The Workflow

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage all \
  --runner local \
  --submit
```

## Check Results

Start with:

- `work/<sample>/summary/sample_summary.tsv`
- `work/<sample>/summary/per_spot_summary.tsv`
- `work/<sample>/qc/saturation/saturation_summary.tsv`

More output guidance: [Results guide](results.md)

## Common Mistakes

- using `scripts/make_cmd.py` instead of `scripts-emseq/make_cmd.py`
- missing `split_barcodes` (the CB whitelist for `call`)
- missing `call_reference_file` or `call_jobs`
- skipping the dry-run
- running host `call` without an indexed `pooled/pooled.host.sorted.bam`

Developer details such as chemistry-specific behavior and stage contracts are documented separately.
