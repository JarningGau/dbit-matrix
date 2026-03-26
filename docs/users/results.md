# Results Guide

This page explains what users should check first after a run.

## Common Checks For TAPS And EMSeq

Start with these files under `work/<sample>/`:

- `summary/sample_summary.tsv`
- `summary/per_spot_summary.tsv`
- `summary/*.heatmap.png`
- `qc/saturation/saturation_summary.tsv`
- `qc/saturation/saturation_curve.png`

Useful methylation outputs:

- `coverage/host/**/*.CG.cov`
- `coverage/host_mito.CG.cov`
- `coverage/<spike_name>.CG.cov` when spike-ins are configured

## Common Checks For RNA

Start with:

- `demux/<sample>.stats.json`
- `solo/`
- `solo/star.Solo.out/Gene/raw/barcodes_pos.tsv`
- `solo/star.Solo.out/Gene/filtered/barcodes_pos.tsv`

## If Something Looks Wrong

Check these first:

- dry-run completed without missing parameters
- input FASTQ paths are correct
- whitelist and reference paths are correct
- the workflow entry script matches the assay type

## Where To Look Next

- stage contracts: [docs/developers/contracts.md](../developers/contracts.md)
- config field definitions: [docs/developers/config-reference.md](../developers/config-reference.md)
- maintainer regression procedures: [docs/maintenance/](../maintenance/)
