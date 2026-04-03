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
- `methscan_*` stages are optional, explicit, and not part of `--stage all` (see [Methscan optional stages](../developers/contracts.md#methscan-optional-stages-methscan_) in contracts)
- no path resolution errors appear
- Slurm dry-run resolves tools from the active `pixi` environment unless explicit `--*-bin` overrides are used

## Optional methscan checks

Run these when changing methscan wiring (`make_cmd`, `scripts/methscan_run.py`, workflow JSON methscan fields, or `envs/methscan`). EMSeq marks this path experimental in code; I/O and layout are still defined in contracts.

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage methscan_prepare \
  --runner local \
  --dry-run

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage methscan_prepare \
  --runner slurm \
  --dry-run

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage methscan_all \
  --runner local \
  --dry-run

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage methscan_all \
  --runner slurm \
  --dry-run
```

Expected results:

- commands expand without missing parameters (`methscan_all` requires `methscan_tss_bed` in the workflow JSON; the test config already sets it)
- Slurm job names use the `emseq_methscan_*` prefix

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

Methscan stage behavior and paths under `work/<sample>/coverage/` are specified in the same contracts (Methscan optional stages).

For command-generation behavior changes, also check [docs/developers/runner-reference.md](../developers/runner-reference.md).
