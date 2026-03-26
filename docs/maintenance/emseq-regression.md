# EMSeq Regression

This document is for maintainers, not first-time users.

## Minimum Checks

Run these before merging an EMSeq workflow change:

```bash
pixi run python scripts-emseq/make_cmd.py --help

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage all \
  --runner local \
  --dry-run

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage all \
  --runner slurm \
  --dry-run
```

## Expected Results

- commands expand without missing parameters
- stage order remains `fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> saturation -> summary`
- `--stage all` covers the fixed mainline only
- no path resolution errors appear
- Slurm dry-run resolves tools from the active `pixi` environment unless explicit `--*-bin` overrides are used

## Optional Real Generation Check

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage all \
  --runner slurm
```

Expected result:

- `work/<sample>/commands/run.sbatch` is generated

## When To Go Deeper

If a stage implementation changed, validate that stage directly and compare outputs against the current contract in [docs/developers/contracts.md](../developers/contracts.md).

For command-generation behavior changes, also check [docs/developers/runner-reference.md](../developers/runner-reference.md).
