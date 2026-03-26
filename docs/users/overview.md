# User Overview

This documentation set is for people who want to run the workflows, not maintain the codebase.

## What This Repo Does

DBiT-Matrix processes DBiT-style spatial omics data and organizes results under `work/<sample>/`.

Supported workflows:

- `TAPS`
- `EMSeq`
- `RNA`

## How To Use The Docs

1. Pick the workflow that matches your assay.
2. Copy the matching example file from [`workflow/`](../../workflow/).
3. Run a dry-run with the workflow entrypoint.
4. Run the workflow.
5. Check outputs in `work/<sample>/`.

## Workflow Guides

- [Environment and inputs](setup.md)
- [TAPS user guide](taps.md)
- [EMSeq user guide](emseq.md)
- [RNA user guide](rna.md)

## Before You Start

You need:

- `pixi install`
- paired-end FASTQ files
- the correct whitelist files for your assay
- reference files required by your workflow

Use the workflow-specific guide for the exact minimum inputs.

## What New Users Should Ignore

You do not need to understand:

- per-stage Slurm job wiring
- internal chunk naming rules
- stage-by-stage contracts
- maintainer regression procedures

Those details live in the developer and maintenance docs.
