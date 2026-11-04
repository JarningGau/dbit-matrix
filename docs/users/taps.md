# TAPS User Guide

Use this guide for the TAPS and TAPS v2 workflows.

## What It Runs

The mainline stages are `fastp_split` through `summary`. Full order and optional stages: [stage contracts](../developers/contracts.md).

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

## TAPS v2

TAPS v2 uses the same mainline and output layout as TAPS. The difference is demux: it matches a C→N-masked methylated-linker and writes conversion QC in `stats.json`.

Start from:

- config templates: [`workflow/dbit_taps_v2_test.json`](../../workflow/dbit_taps_v2_test.json), [`workflow/dbit_taps_v2_slurm.json`](../../workflow/dbit_taps_v2_slurm.json)
- entry script: `scripts-v2/make_cmd.py`

```bash
pixi run python scripts-v2/make_cmd.py \
  --workflow-config workflow/dbit_taps_v2_test.json \
  --stage all \
  --runner local \
  --dry-run
```

Optional demux filter: `require_c_all_t` (empty = off). Details: [TAPS v2 demux contract](../developers/contracts.md#taps-v2-demux-scriptsv2extract_bcpy) and [config reference](../developers/config-reference.md#taps-fields).

## Optional: 5hmC Oligo QC (`hmc_spike`)

TAPS and TAPS v2 can score a synthetic 5hmC oligo after `saturation` and before `summary`. EMSeq does not support this stage.

Set `hmc_spike.enabled` to `true` in the workflow JSON (see [`workflow/dbit_taps_test.json`](../../workflow/dbit_taps_test.json)). Outputs land under `qc/hmc_spike/`.

Details: [hmc_spike contract](../developers/contracts.md#hmc_spike) and [config reference](../developers/config-reference.md#taps-fields).

## Check Results

Start with:

- `work/<sample>/summary/sample_summary.tsv`
- `work/<sample>/summary/per_spot_summary.tsv`
- `work/<sample>/qc/saturation/saturation_summary.tsv`

More output guidance: [Results guide](results.md)

## Common Mistakes

- using the wrong entry script (`scripts/make_cmd.py` vs `scripts-v2/make_cmd.py`)
- missing whitelist paths
- missing reference paths
- skipping the dry-run
- running host `call` without an indexed `pooled/pooled.host.sorted.bam`

Developer details such as stage contracts and full config fields are documented separately.
