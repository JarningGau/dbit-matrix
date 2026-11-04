# AGENTS.md

Use these rules when working in this repository.

## Core Principles

- Start from first principles.
- Prefer the simplest solution that preserves the required contract.
- Build the MVP first. Extend only when the MVP is stable.
- Make one main change per iteration.

## Current Workflow Status

- TAPS, TAPS v2, and EMSeq share this mainline:
  `fastp_split -> demux_extract_bc -> align -> pool -> mbias -> call -> saturation -> summary`
- TAPS / TAPS v2 may insert `hmc_spike` after `saturation` and before `summary` when `hmc_spike.enabled` is true. It depends only on demux `*.spike-in.fq.gz` (not a prerequisite for `align`). EMSeq does not implement `hmc_spike`.
- `--stage all` generates `work/<sample>/commands/run.sh` or `run.sbatch` in stage order.
- `split` is not in `--stage all` or `STAGE_CHOICES`. `scripts/split_bams.py` exists but is not a workflow stage.
- `aggregate` and methscan optional stages (`methscan_prepare`, `methscan_filter`, `methscan_profile`, `methscan_smooth`, `methscan_scan`, `methscan_matrix`, `methscan_all`) are optional explicit stages and are not part of `--stage all`.
- `scripts-emseq/mbias.py` is the EMSeq `mbias` implementation and is different from TAPS `scripts/mbias.py`.
- EMSeq `call` is pooled pysam (`scripts-emseq/methy_caller.py` / `scripts-emseq/call.py`), not biscuit pileup.

## Engineering Rules

- Keep single-stage scripts thin, explicit, and predictable.
- Put workflow parameters in `workflow/*.json` whenever possible.
- Every stage must support `--dry-run`.
- Update `pixi.lock` whenever dependencies change.
- Do not change an input/output contract silently.

## Release version

When cutting a release, update these in the same change (single source of truth for CLI: `scripts/_version.py`):

| Artifact | Action |
|----------|--------|
| `scripts/_version.py` | Set `__version__` to the new semver (drives `--version` on all workflow drivers). |
| `README.md` | Update the “Current version” line. |
| `docs/maintenance/changelog.md` | Move `## Current` items into `## X.Y.Z`, reset `## Current` (per existing style). |

Do **not** treat `envs/methscan/pixi.toml` `version` as the repository release; that field is for the methscan sub-environment only.

## Entrypoints

- TAPS workflow driver: `scripts/make_cmd.py`
- TAPS v2 workflow driver: `scripts-v2/make_cmd.py` (methylated-linker demux; reuses most `scripts/*.py`)
- EMSeq workflow driver: `scripts-emseq/make_cmd.py`
- RNA workflow driver: `scripts-rna/make_cmd.py`
- Single-stage implementations live in:
  - `scripts/*.py`
  - `scripts-v2/*.py` (forked demux/summary only)
  - `scripts-emseq/*.py`
  - `scripts-rna/*.py`

## Slurm Rules

- Use the `pixi` environment. Do not rely on `module load`.
- Keep Slurm settings stage-specific. Do not reuse one generic resource block for unrelated stages.
- Do not add unrelated generic options to sbatch templates.
- Required stage layout:
  - `demux_extract_bc`: one sbatch per chunk
  - `align`: one sbatch per chunk; run spike-in before host inside the script
  - `pool`: separate host and spike-in jobs
  - `mbias`: split host and spike jobs; use `slurm.mbias.host` and `slurm.mbias.spike`
  - `call`: split host and spike jobs; use `slurm.call.host` and `slurm.call.spike`
  - `saturation`: one sbatch; use `slurm.saturation`
  - `hmc_spike` (TAPS / TAPS v2, when enabled): one sbatch; use `slurm.hmc_spike`
  - `summary`: one sbatch; use `slurm.summary`

## Key Contracts

- `demux` writes host demux FASTQ and spike-in FASTQ:
  - `*.demux.fq.gz`
  - `*.spike-in.fq.gz`
- `demux` stats must include retention information and reject reasons.
- `align` reads host input from `*.demux.fq.gz` and spike-in input from `*.spike-in.fq.gz`.
- `align` writes:
  - host: `<chunk>.cb.bam`
  - spike-in: `<chunk>.<spike_name>.bam`
- `spike_in_index` must support either a JSON object or a `NAME=INDEX` list.
- `pool` writes:
  - host: `pooled/pooled.host.sorted.bam` (coordinate-sorted, indexed; no `-t CB`)
  - spike-in: `pooled/pooled.<spike_name>.sorted.bam` (mapped reads only)
- Intermediate CB is concatenated `<bcB><bcA>` (no `+`); host cov is `coverage/host/<cb[0:3]>/<cb>.CG.cov`.
- `call` writes `coverage/per_spot_read_counts.tsv` (`barcode`, `reads`). Smoke is `call_smoke` (top 16 CBs by reads).
- `split_barcodes` is the CB whitelist for `call`, not a split-stage input.
- Summary uses `barcode_positions` for `array_row`/`array_col`; summary `barcode` is `<bcB><bcA>`.

## Documentation Rules

- Write new documentation in English.
- Use the `docs/` tree as the source of truth for active docs.
- Update the right doc when behavior changes:
  - `README.md`: homepage, workflow status, top-level navigation
  - `docs/users/`: user-facing run guides
  - `docs/developers/contracts.md`: stage contracts
  - `docs/developers/config-reference.md`: config fields and runner rules
  - `docs/maintenance/`: maintainer regression procedures
- Do not duplicate one normative fact across multiple docs. Link to the source-of-truth page instead.

## Required Checks

Before finishing a change, run at least:

- relevant CLI `--help` and `--version`
- relevant workflow `--dry-run`

If contracts changed, also update the matching docs in the same change.

## Commit Style

- Use [Conventional Commits](https://www.conventionalcommits.org/): `type: short description`
- Common types: `feat`, `fix`, `refactor`, `docs`, `chore`, `data`, `analysis`
- Write commit messages in English
- Write the subject line in imperative mood, lowercase, with no trailing period
- Keep the subject line under 72 characters
- One commit should focus on one reason
- Use the body only when the why is not obvious from the subject
- Separate the body from the subject with a blank line
- If the contract changed, say so explicitly in the body
- Use `BREAKING:` in the body for breaking changes

## Safety

- Do not overwrite or revert user changes without explicit permission.
- If unexpected changes conflict with your task, stop and confirm before proceeding.
