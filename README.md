# DBiT-Matrix

DBiT-Matrix is a workflow repository for DBiT-style spatial omics data processing. It currently provides three entrypoints:

- `TAPS` for methylation workflows
- `EMSeq` for methylation workflows
- `RNA` for transcriptomics workflows

All workflows are configured with [`workflow/*.json`](workflow/) and run through `pixi`.

## Current Status

- `TAPS`: main methylation workflow available
- `EMSeq`: main methylation workflow available
- `RNA`: current scope is `demux_extract_bc -> align`

## Install

```bash
pixi install
```

## Choose A Workflow

| Workflow | Entrypoint | Start here |
| --- | --- | --- |
| `TAPS` | `scripts/make_cmd.py` | [TAPS user guide](docs/users/taps.md) |
| `EMSeq` | `scripts-emseq/make_cmd.py` | [EMSeq user guide](docs/users/emseq.md) |
| `RNA` | `scripts-rna/make_cmd.py` | [RNA user guide](docs/users/rna.md) |

General orientation for users: [User overview](docs/users/overview.md)

## First Run Pattern

Start with a dry-run before submitting a real job:

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage all \
  --runner local \
  --dry-run
```

Use the workflow-specific guide above to select the correct entry script and config template.

## Main Outputs

Outputs are written under `work/<sample>/`. Common downstream files include:

- `summary/sample_summary.tsv`
- `summary/per_spot_summary.tsv`
- `qc/saturation/saturation_summary.tsv`
- workflow-specific files under `coverage/`, `qc/`, `solo/`, and `commands/`

User-facing output checks: [Results guide](docs/users/results.md)

## Documentation

### Users

- [User overview](docs/users/overview.md)
- [TAPS user guide](docs/users/taps.md)
- [EMSeq user guide](docs/users/emseq.md)
- [RNA user guide](docs/users/rna.md)
- [Results guide](docs/users/results.md)

### Developers

- [Architecture](docs/developers/architecture.md)
- [Stage contracts](docs/developers/contracts.md)
- [Config reference](docs/developers/config-reference.md)
- [Doc system](docs/developers/doc-system.md)

### Maintenance

- [TAPS regression](docs/maintenance/taps-regression.md)
- [EMSeq regression](docs/maintenance/emseq-regression.md)
- [RNA regression](docs/maintenance/rna-regression.md)

### Legacy

Legacy Chinese docs are preserved in [`doc/`](doc/). New work should target [`docs/`](docs/).
