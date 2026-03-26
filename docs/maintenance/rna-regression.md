# RNA Regression

This document is for maintainers, not first-time users.

## Minimum Checks

Run these before merging an RNA workflow change:

```bash
pixi run python scripts-rna/make_cmd.py --help
pixi run python scripts-rna/extract_bc.py --help
pixi run python scripts-rna/align.py --help

pixi run python scripts-rna/make_cmd.py \
  --workflow-config workflow/dbit_rna_test.json \
  --stage all \
  --runner local \
  --dry-run

pixi run python scripts-rna/make_cmd.py \
  --workflow-config workflow/dbit_rna_test.json \
  --stage all \
  --runner slurm \
  --dry-run
```

## Expected Results

- commands expand without missing parameters
- stage order remains `demux_extract_bc -> align`
- no path resolution errors appear
- outputs point to `demux/`, `solo/`, and `commands/`

## Optional Real Generation Check

```bash
pixi run python scripts-rna/make_cmd.py \
  --workflow-config workflow/dbit_rna_test.json \
  --stage all \
  --runner local
```

Expected result:

- the generated command set targets `demux/` and `solo/`

## When To Go Deeper

If a stage implementation changed, validate that stage directly and compare outputs against the current contract in [docs/developers/contracts.md](../developers/contracts.md).

For command-generation behavior changes, also check [docs/developers/runner-reference.md](../developers/runner-reference.md).
