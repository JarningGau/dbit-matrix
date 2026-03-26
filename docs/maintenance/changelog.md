# Changelog

This page is the active release history for the repository.

## 1.6.0

- added the standalone RNA workflow with:
  `demux_extract_bc -> align`
- added RNA stage scripts and workflow generation for local and Slurm runs
- pinned RNA Slurm logs under `work/<sample>/logs/`
- added the experimental `aggregate` stage for TAPS and EMSeq host coverage outputs
- added the experimental `methscan_prepare` stage for host coverage outputs
- exposed `methscan_prepare` `chunksize` through workflow config
- wired EMSeq calling trimming to biscuit pileup `-5/-3` workflow settings
- removed the deprecated `scripts-emseq/extract_bc.py` wrapper and unified EMSeq demux behavior on `scripts/extract_bc.py`
- reduced memory pressure during CG aggregation and allowed reruns to reuse existing aggregated outputs
- added `cg_span` histogram generation
- rebuilt the documentation system in English under `docs/`
- turned `README.md` into the main homepage and navigation page
- moved maintainer guidance behind dedicated maintenance pages
- simplified `AGENTS.md` and aligned commit rules with Conventional Commits

## 1.5.0

- added the standalone EMSeq workflow entrypoint
- completed the EMSeq mainline:
  `fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> saturation -> summary`
- added EMSeq stage scripts and workflow examples
- added `--stage all` support for ordered local and Slurm command generation

## 1.4.0

- added the `saturation` stage
- added `saturation_rate` to `summary/sample_summary.tsv`
- added workflow config and Slurm resource support for `saturation`

## 1.3.0

- changed `call` to stream `.CG.cov` output by batch
- aligned trimming semantics between `call` and `mbias`
- reworked summary reads metrics and heatmap inputs
- added a Slurm workflow profile

## 1.2.0

- changed `all + slurm + submit` to client-side DAG submission
- pre-expanded `demux_extract_bc` and `align` Slurm jobs from split-part counts
- clarified spike-in ownership in `mbias` and `call`

## 1.1.3

- clarified user-facing `call` outputs in the docs
- fixed the Slurm `split` stage to keep `split_bams` and `sort` separate
- updated TAPS test workflow Slurm defaults

## 1.1.2

- preferred tool resolution from the active `pixi` environment
- removed default `module load fastp` from Slurm `fastp_split`
- added `05_split_submit.sh` as a split-stage submission helper

## 1.1.1

- slowed the default progress print interval in `scripts/extract_bc.py`
- updated the TAPS test workflow to generate both host and spike-in `mbias` and `call`

## 1.1.0

- added three summary heatmap outputs
- added summary documentation and testing notes

## 1.0.0

- established the original mainline:
  `fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> summary`
- supported both local and Slurm execution with `pixi`
- delivered the first closed user-facing summary and `mbias` outputs
