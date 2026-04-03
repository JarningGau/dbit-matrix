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
- `aggregate` and methscan optional stages (`methscan_prepare`, `methscan_filter`, `methscan_profile`, `methscan_smooth`, `methscan_scan`, `methscan_matrix`, `methscan_all`) are explicit and are not part of `all`

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

### Methscan optional stages

- available for TAPS and EMSeq only
- generated commands invoke `scripts/methscan_run.py` (see `contracts.md`) with `--work-path work/<sample>`
- `methscan_prepare` reads `coverage/host/**/*.CG.cov` and writes `coverage/host_prepare/`
- `methscan_filter` reads `coverage/host_prepare/` and writes `coverage/filter/`
- `methscan_profile` requires `methscan_tss_bed`; default profile input dir is `coverage/host_prepare/`
- `methscan_smooth` reads `coverage/filter/` (writes `coverage/filter/smoothed/` per methscan)
- `methscan_scan` reads `coverage/filter/` and writes VMRs (default `coverage/VMRs.bed`)
- `methscan_matrix` reads VMRs bed + `coverage/filter/` (default matrix dir `coverage/matrix_VMR/`)
- `methscan_all` runs the six steps in one script; requires `methscan_tss_bed`
- uses the `envs/methscan` pixi workspace unless `methscan_pixi_manifest` overrides it
