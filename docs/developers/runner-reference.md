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

- `TAPS`, `TAPS v2`, and `EMSeq` expand `all` in stage order from [contracts](contracts.md); `hmc_spike` is included only when `hmc_spike.enabled` is true (TAPS / TAPS v2)
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
- stage dependencies use `afterok` edges that follow data deps (not only adjacent `STAGE_SEQUENCE` order)
- when `hmc_spike` is enabled, its jobs depend on `demux_extract_bc`; `summary` depends on both `saturation` and `hmc_spike`
- nested `sbatch` on compute nodes is not part of the model

In Slurm mode, `mbias` and `call` remain separate host and spike jobs. See [config reference](config-reference.md#slurm-structure).

## Advanced Stages

### `aggregate`

- available for TAPS and EMSeq only
- reads `coverage/host/**/*.CG.cov`
- writes `coverage/aggregated_cg_by_id.tsv` and `coverage/aggregated_cg_by_pos.tsv`

### Methscan optional stages

- available for TAPS and EMSeq only
- generated commands invoke `scripts/methscan_run.py` (see `contracts.md`) with `--work-path work/<sample>`
- `methscan_prepare` reads `coverage/host/**/*.CG.cov` and writes `methscan/compact/`
- `methscan_filter` reads `methscan/compact/` and writes `methscan/filter/`
- `methscan_profile` requires `methscan_tss_bed`; default profile input dir is `methscan/compact/`
- `methscan_smooth` reads `methscan/filter/` (writes `methscan/filter/smoothed/` per methscan)
- `methscan_scan` reads `methscan/filter/` and writes VMRs (default `methscan/VMRs.bed`)
- `methscan_matrix` reads VMRs bed + `methscan/filter/` (default matrix dir `methscan/matrix/`)
- `methscan_all` runs the six steps in one script; requires `methscan_tss_bed`
- uses the `envs/methscan` pixi workspace unless `methscan_pixi_manifest` overrides it
