# AGENTS.md

Use these rules when working in this repository.

## Core Principles

- Start from first principles.
- Prefer the simplest solution that preserves the required contract.
- Build the MVP first. Extend only when the MVP is stable.
- Make one main change per iteration.

## Current Workflow Status

- EMSeq mainline is fixed:
  `fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> saturation -> summary`
- EMSeq `--stage all` generates `work/<sample>/commands/run.sh` or `run.sbatch` in stage order.
- `aggregate` and methscan optional stages (`methscan_prepare`, `methscan_filter`, `methscan_profile`, `methscan_smooth`, `methscan_scan`, `methscan_matrix`, `methscan_all`) are optional explicit stages and are not part of `--stage all`.
- `scripts-emseq/mbias.py` is the EMSeq `mbias` implementation and is different from TAPS `scripts/mbias.py`.

## Engineering Rules

- Keep single-stage scripts thin, explicit, and predictable.
- Put workflow parameters in `workflow/*.json` whenever possible.
- Every stage must support `--dry-run`.
- Update `pixi.lock` whenever dependencies change.
- Do not change an input/output contract silently.

## Entrypoints

- TAPS workflow driver: `scripts/make_cmd.py`
- EMSeq workflow driver: `scripts-emseq/make_cmd.py`
- RNA workflow driver: `scripts-rna/make_cmd.py`
- Single-stage implementations live in:
  - `scripts/*.py`
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
  - `split`: `split_bams` and `sort` must be separate sbatch jobs with dependencies
  - `mbias`: split host and spike jobs; use `slurm.mbias.host` and `slurm.mbias.spike`
  - `call`: split host and spike jobs; use `slurm.call.host` and `slurm.call.spike`
  - `saturation`: one sbatch; use `slurm.saturation`
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
  - host: `pooled/pooled.byCB.bam`
  - spike-in: `pooled/pooled.<spike_name>.sorted.bam`
- `split` reads `pooled/pooled.byCB.bam` and parses spots from `CB:Z:<x>+<y>`.
- `split --smoke` writes at most 16 non-empty spot BAMs.

## Documentation Rules

- Write new documentation in English.
- Use the new `docs/` tree as the source of truth for active docs.
- Treat `doc/` as legacy reference unless a task explicitly targets it.
- Update the right doc when behavior changes:
  - `README.md`: homepage, workflow status, top-level navigation
  - `docs/users/`: user-facing run guides
  - `docs/developers/contracts.md`: stage contracts
  - `docs/developers/config-reference.md`: config fields and runner rules
  - `docs/maintenance/`: maintainer regression procedures
- Do not duplicate one normative fact across multiple docs. Link to the source-of-truth page instead.

## Required Checks

Before finishing a change, run at least:

- relevant CLI `--help`
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
