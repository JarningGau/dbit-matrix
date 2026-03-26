# Runner Reference

This page documents command-generation behavior that is too detailed for user guides.

## Common Command Pattern

All workflow drivers follow this shape:

```bash
pixi run python <entry-script> \
  --workflow-config <workflow.json> \
  --stage <stage-or-all> \
  --runner <local-or-slurm> \
  --dry-run
```

## `all`

- `TAPS` and `EMSeq` expand `all` to:
  `fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> saturation -> summary`
- `RNA` expands `all` to:
  `demux_extract_bc -> align`
- `aggregate` and `methscan_prepare` are explicit optional stages and are not part of `all`

## `local`

- writes per-stage `.sh` files under `work/<sample>/commands/`
- writes `run.sh` for `--stage all`
- `--submit` executes the generated shell entrypoint

## `slurm`

- writes per-stage or per-chunk `.sbatch` files under `work/<sample>/commands/`
- writes `run.sbatch` for `--stage all`
- `--submit` submits the dependency graph from the client side
- stage dependencies use `afterok`
- nested `sbatch` on compute nodes is not part of the model

## Required Split Layout

In Slurm mode, `split` is always two jobs:

- `split_bams`
- `sort`

The sort job must depend on the split job.

## Advanced Stages

### `aggregate`

- available for TAPS and EMSeq only
- reads `coverage/host/**/*.CG.cov`
- writes `coverage/aggregated_cg_by_id.tsv` and `coverage/aggregated_cg_by_pos.tsv`

### `methscan_prepare`

- available for TAPS and EMSeq only
- reads `coverage/host/**/*.CG.cov`
- writes `coverage/host_prepare/`
- uses the `envs/methscan` pixi workspace unless overridden
