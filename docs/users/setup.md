# Environment And Inputs

This page is for first-time users who want a quick checklist before running a workflow.

## Install

```bash
pixi install
```

Optional environment checks:

```bash
pixi run which python
pixi run python --version
```

## What You Need

All workflows need:

- paired-end FASTQ files
- the correct barcode whitelist files
- a workflow JSON copied from `workflow/`

Methylation workflows also need:

- host reference files for alignment and calling
- optional spike-in references when used

RNA also needs:

- `star_genome_dir`
- `gtf`
- `umi_left`
- `umi_len`

## Starter Configs

Use one of these as a starting point:

- [`workflow/dbit_taps_test.json`](../../workflow/dbit_taps_test.json)
- [`workflow/dbit_emseq_test.json`](../../workflow/dbit_emseq_test.json)
- [`workflow/dbit_rna_test.json`](../../workflow/dbit_rna_test.json)

## Minimal Library Structure Assumption

For TAPS and EMSeq demux, the workflow only relies on:

- `linker_bc`
- `insert_left`

The pipeline does not require a full fixed template beyond those anchors.

## Next Step

After your environment and inputs are ready, continue with:

- [TAPS user guide](taps.md)
- [EMSeq user guide](emseq.md)
- [RNA user guide](rna.md)
