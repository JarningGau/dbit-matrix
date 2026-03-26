# TAPS User Guide

Use this guide for the TAPS workflow.

## What It Runs

The main pipeline is:

`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> saturation -> summary`

## Before You Start

Prepare:

- `R1` and `R2` FASTQ files
- `barcode1_whitelist`
- `barcode2_whitelist`
- host reference files for alignment and calling
- optional spike-in reference files

Start from:

- config template: [`workflow/dbit_taps_test.json`](../../workflow/dbit_taps_test.json)
- entry script: `scripts/make_cmd.py`

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
- `bwa_index`
- `call_reference_file`

## Dry-Run

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner local \
  --dry-run
```

If you use Slurm, change `--runner local` to `--runner slurm`.

## Run The Workflow

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
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

- using the wrong entry script
- missing whitelist paths
- missing reference paths
- skipping the dry-run

Developer details such as stage contracts and full config fields are documented separately.
