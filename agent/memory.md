# Workflow Development Memory

## Part I. Design Preferences

### 1. Pipeline Architecture

- DO define the mainline pipeline in one explicit ordered stage list.
- DO implement `--stage all` by iterating that ordered stage list.
- DO keep optional analyses as explicit opt-in stages instead of silently including them in `--stage all`.
- DO keep stage execution order derived from one canonical source instead of duplicating the order in multiple code paths.

### 2. HPC / Slurm

- DO keep single-stage scripts thin and execution-only.
- DO generate runner-specific script files for both `local` and `slurm` execution.
- DO keep Slurm resource settings stage-specific instead of forcing unrelated stages to share one generic resource block.
- DO allow step-level Slurm overrides when one stage expands into multiple jobs.

### 3. Environment & Tooling

- DO use `pixi` as the default pipeline environment entrypoint.
- DO prefer executables from the current Python environment before falling back to bare command names.
- DO quote generated shell command arguments with shell-safe escaping.
- DO generate shell scripts with `set -euo pipefail`.

## Part II. Hard Constraints

### 1. Pipeline Behavior

- DO support both `local` and `slurm` runners through the same workflow driver.
- DO validate stage-required settings before generating commands for that stage.
- DO implement `--dry-run` so it prints the generated command and script path without writing files or submitting jobs.
- DO NOT write workflow scripts when `--dry-run` is set.
- DO NOT submit jobs when `--dry-run` is set.

### 2. HPC / Slurm

- DO submit upstream-to-downstream Slurm dependencies with `sbatch --dependency=afterok:<jobid>`.
- DO keep Slurm orchestration in a client-side driver or submit helper rather than inside submitted stage scripts.
- DO NOT call `sbatch` from inside a stage script that is itself submitted by `sbatch`.
- DO NOT rely on `module load` for pipeline execution.

### 3. Data Handling

- DO validate required input paths before writing stage scripts.
- DO validate required upstream workdir outputs before generating a downstream stage, unless an explicit orchestration flag skips those checks.
- DO normalize config fields that accept multiple input shapes into one internal representation before building commands. Example: normalize `spike_in_index` from either a JSON object or a `NAME=INDEX` list into a `NAME=INDEX` list.
- DO reject malformed normalized config entries. Example: reject `spike_in_index` values without both a non-empty name and a non-empty index.
- DO NOT defer obvious missing-input failures until after subprocess fan-out or job submission.
- DO NOT make `conda` or shell modules the default execution path when a `pixi` workflow path exists.

## Part III. Anti-patterns

- DO reject nested Slurm submission patterns because they cause job explosion. Example: if `01_align.sbatch` is submitted by `sbatch`, `01_align.sbatch` must not run another `sbatch`.
- DO reject workflow designs that silently bundle optional analyses into `--stage all`.
- DO reject workflow designs that skip preflight validation and rely on downstream subprocesses or cluster jobs to discover missing required inputs.
