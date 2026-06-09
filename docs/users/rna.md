# RNA User Guide

Use this guide for the RNA workflow.

## Library Structure

R1: BC2 (8 bp) + linker2 + BC1 (8 bp) + linker1 + UMI (10 bp)

R2: cDNA

linker1 (30 bp): GTGGCCGATGTTTCGCATCGGCGTACGACT (umi_left)

linker2 (30 bp): ATCCACGTGCTTGAGAGGCCAGAGCATTCG (linker_bc)

## What It Runs

The current RNA scope is:

`demux_extract_bc -> align`

The RNA entrypoint does not use chunked processing.

## Before You Start

Prepare:

- `R1` and `R2` FASTQ files
- `barcode1_whitelist`
- `barcode2_whitelist`
- `linker_bc`
- `umi_left`
- `umi_len`
- `star_genome_dir`
- `gtf`

Start from:

- config template: [`workflow/dbit_rna_test.json`](../../workflow/dbit_rna_test.json)
- entry script: `scripts-rna/make_cmd.py`

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
- `umi_left`
- `umi_len`
- `star_genome_dir`
- `gtf`
- `solo_features` (must include `Gene` for `barcodes_pos.tsv` generation; default: `Gene GeneFull`)

## Dry-Run

```bash
pixi run python scripts-rna/make_cmd.py \
  --workflow-config workflow/dbit_rna_test.json \
  --stage all \
  --runner local \
  --dry-run
```

If you use Slurm, change `--runner local` to `--runner slurm`.

## Run The Workflow

```bash
pixi run python scripts-rna/make_cmd.py \
  --workflow-config workflow/dbit_rna_test.json \
  --stage all \
  --runner local \
  --submit
```

## Check Results

Start with:

- `work/<sample>/demux/<sample>.stats.json`
- `work/<sample>/solo/`
- `work/<sample>/solo/star.Solo.out/Gene/raw/barcodes_pos.tsv`

## Common Mistakes

- using a methylation entry script instead of `scripts-rna/make_cmd.py`
- missing `star_genome_dir` or `gtf`
- missing `umi_left` or `umi_len`
- expecting methylation coverage outputs from the RNA workflow
