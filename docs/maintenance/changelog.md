# Changelog

This page is the active release history for the repository.

## Current

- **BREAKING saturation auto knee:** auto HQ cutoff also excludes the bottom 10% of barcode ranks from the slope search (still excludes top 5%, then leftmost half-steep). Stops a steeper right-hand tail from dominating `min_slope`. Re-run `saturation` then `summary` when using `auto`. See [contracts](../developers/contracts.md#saturation).
- **Summary call sidecars (fast path):** host `call` writes `coverage/host_bam_funnel.json` (same BAM pass as `per_spot_read_counts.tsv`) and `coverage/per_spot_meth_stats.tsv` (while shuffling host cov). `summary` prefers both for host funnel + per-spot methylation; missing/invalid sidecars warn and fall back to the previous slow host BAM/cov scan. Re-run host `call` then `summary` for the minutes-level path. See [contracts](../developers/contracts.md#call) and [contracts](../developers/contracts.md#summary).
- **BREAKING hmc_spike motif match:** locate requires exact non-`C` bases and `[CT]` at every `C` in `sequence` (no backbone `edit_distance`). Re-run `hmc_spike` then `summary`. See [contracts](../developers/contracts.md#hmc_spike).
- **BREAKING summary `is_HQ_spot` and column order:** `per_spot_summary.tsv` gains `is_HQ_spot` (`1`/`0`/`NA` from saturation `reads_threshold`) and reorders columns to identity → `reads`/`is_HQ_spot` → methylation → site counts. `sample_summary.tsv` reorders to `sample_id`/`saturation_rate` → read funnel → (v2 linker QC) → host/mito/spike methylation → optional `*_hmc_spike_*`. Re-run `summary`. See [contracts](../developers/contracts.md#summary).
- **BREAKING host_mito BAM source:** `call` builds `coverage/host_mito.*.cov` from full indexed `pooled/pooled.host.sorted.bam` (mito contigs only), not from `qc/mbias/host.subsampled.sorted.bam`. Re-run `call` then `summary`. Removed EMSeq `call_host_subsample_fraction` / `call_host_subsample_seed` and call-stage subsample/samtools CLI. See [contracts](../developers/contracts.md#call).
- **mbias host subsample unbiased:** [`scripts/host_subsample_bam.py`](../../scripts/host_subsample_bam.py) drops the coordinate-ordered `head -n` cap. Size is limited by lowering the `samtools view -s` fraction from `idxstats` mapped counts (`HOST_SUBSAMPLE_MAX_READS`). Re-run `mbias` to refresh `qc/mbias/host.subsampled.bam`. See [contracts](../developers/contracts.md#mbias).
- **BREAKING mbias host subsample artifact:** mbias no longer writes `qc/mbias/host.subsampled.sorted.bam` / `.bai`; host M-bias reads unsorted `qc/mbias/host.subsampled.bam`. See [contracts](../developers/contracts.md#mbias).

## 4.0.0

- **BREAKING CB / coverage / summary identity:** Intermediate CB is `<bcB><bcA>` (no `+`); host cov paths and summary columns use `barcode` / `array_row` / `array_col`. Re-run `call` then `summary`. See [contracts](../developers/contracts.md#call) and [contracts](../developers/contracts.md#summary).
- **BREAKING spike pooled BAM:** `pooled/pooled.<spike_name>.sorted.bam` keeps mapped reads only. Re-run `pool`. See [contracts](../developers/contracts.md#pool).
- **BREAKING summary barcode / 5mCG / 5mCH / float precision:** `barcode` drops the `host` suffix; `mean_methylation` becomes `mean_5mCG` / `mean_5mCH`; non-integer summary floats use two decimals. Re-run `summary`. See [contracts](../developers/contracts.md#summary).
- **BREAKING saturation HQ threshold:** default is `"auto"` (was `1e6`); auto knee excludes the top 5% of ranks and picks the leftmost half-steep slope. Writes `reads_knee.png`. Re-run `saturation` then `summary`. See [contracts](../developers/contracts.md#saturation).
- **BREAKING hmc_spike rates:** hmC→T rates are percent-scale (0–100). Re-run `hmc_spike` then `summary`. See [contracts](../developers/contracts.md#hmc_spike).
- **TAPS / TAPS v2 `hmc_spike`:** optional 5hmC oligo QC after `saturation` when `hmc_spike.enabled` is true. See [contracts](../developers/contracts.md#hmc_spike).
- **Summary CA/CC/CT and HQ heatmaps:** adds `mean_5mCA` / `mean_5mCC` / `mean_5mCT` fields and heatmaps; methylation heatmaps plot HQ spots only; origin is top-right; Nature-style figure size. See [contracts](../developers/contracts.md#summary).
- **Docs:** user guides cover TAPS v2 and `hmc_spike`; stage order is linked from contracts instead of copied; `docs/superpowers` removed.

## 3.0.0

- **BREAKING host pooled BAM:** Host pooled BAM is `pooled/pooled.host.sorted.bam`, coordinate-sorted and indexed. `pooled.byCB.bam` is no longer written. See [contracts](../developers/contracts.md#pool).
- **BREAKING splitless mainline:** `split` is removed from `--stage all` and from `STAGE_CHOICES`. Mainline is `fastp_split -> demux_extract_bc -> align -> pool -> mbias -> call -> saturation -> summary`. See [contracts](../developers/contracts.md).
- **BREAKING per-spot read counts:** `saturation` / `summary` read `coverage/per_spot_read_counts.tsv`, written by `call`.
- **BREAKING call_smoke:** `split_smoke` is replaced by `call_smoke` (deterministic top 16). See [config reference](../developers/config-reference.md).
- **BREAKING EMSeq call:** EMSeq host/mito/spike calling no longer uses biscuit pileup or VCF intermediates. EMSeq CG percentages and site sets will move. Caller is pooled pysam (`scripts-emseq/methy_caller.py` / `scripts-emseq/call.py`). Compare against the frozen biscuit baseline with `scripts-emseq/compare_cov_to_biscuit.py` and `work/test-DNAme-EMSeq` (do not overwrite).
- **EMSeq biscuit compare (16-spot roster on sibling work tree `work/test-DNAme-EMSeq-pooled-call`):** shared_sites=1160, only_baseline=885, only_new=42, covered_shared=1160, pearson_meth=0.921, pearson_unmeth=0.957, mean_abs_pct_delta=0.73, disagree_majority_fraction=0.0069 (pass bar `< 0.5`).
- **BREAKING EMSeq trim keys:** EMSeq trim keys `call_left_trimming` / `call_right_trimming` are removed; use `call_trimming` like TAPS.
- **BREAKING EMSeq `.CG.cov` coordinates:** EMSeq `.CG.cov` positions are pysam 0-based, not biscuit's transformed coordinate. `pct` has two decimal places.
- **BREAKING TAPS/TAPS v2 call:** host methylation calling reads
  `pooled/pooled.host.sorted.bam` with interval-parallel pileup instead of
  per-spot split BAMs. CG/CH share one pass. CH neighbor/strand rules
  match the template TAPS caller. Defaults: MAPQ 10, `max_depth` 50000.
  Removed `scripts/methy_caller_CH.py` and `call_ch_caller_script`.
  Per-spot `.cov` paths are unchanged. See [contracts](../developers/contracts.md#call).
- **Summary host valid-read mito/nuclear split:** `sample_summary.tsv` gains `host_valid_reads_mito` and `host_valid_reads_nuclear`, partitioning existing `host_valid_reads` from `pooled/pooled.host.sorted.bam` by mitochondrial contig membership (`--mito-chromosomes`, default `chrM`; wired from `call_mito_chromosomes` in TAPS/EMSeq/v2 `make_cmd`). See [contracts](../developers/contracts.md#summary).
- **Host mbias nuclear chromosomes only:** host M-bias CG counting in [`scripts/mbias.py`](../../scripts/mbias.py) and [`scripts-emseq/mbias.py`](../../scripts-emseq/mbias.py) now requires `--chromosomes` and skips contigs outside that whitelist, so `chrM` no longer dilutes `host.mbias.tsv`. TAPS/EMSeq/v2 `make_cmd` pass existing `call_chromosomes` into host mbias. See [contracts](../developers/contracts.md#mbias).
- **Saturation plot/model update:** [`scripts/saturation.py`](../../scripts/saturation.py) now plots median unique CpGs per spot with IQR error bars, uses a `Sequencing depth (Gbp)` x-axis derived from `shard_fastq/fastp.json` (auto-detected; falls back to raw coverage fraction when the report is absent), and selects between a through-origin linear model and the exponential saturation curve via `--linear-r2-threshold` (default `0.99`). The saturation summary TSV gains an `extrapolation_model` column. HQ spot filtering and the unique-CpG metric are unchanged. See [contracts](../developers/contracts.md#saturation).
- **TAPS per-target call trimming:** the `call` stage now supports independent r1/r2 left/right trimming for host and spike-in via a nested `call_trimming: {host: {...}, spike: {...}}` config in both [`scripts/make_cmd.py`](../../scripts/make_cmd.py) and [`scripts-v2/make_cmd.py`](../../scripts-v2/make_cmd.py). Missing fields fall back to the legacy flat `call_r*_trimming` keys, preserving existing behavior. On the local runner, `call_mode: all` with differing host/spike trims emits sequential `call.py` invocations; the Slurm runner already splits host / per-spike sbatch jobs, which now carry per-target trim flags. Example configs updated: [`workflow/dbit_taps_test.json`](../../workflow/dbit_taps_test.json), [`workflow/dbit_taps_slurm.json`](../../workflow/dbit_taps_slurm.json), [`workflow/dbit_taps_v2_test.json`](../../workflow/dbit_taps_v2_test.json), [`workflow/dbit_taps_v2_slurm.json`](../../workflow/dbit_taps_v2_slurm.json). See [config reference](../developers/config-reference.md#taps-fields).
- **TAPS v2 methylated-linker demux:** added parallel driver [`scripts-v2/make_cmd.py`](../../scripts-v2/make_cmd.py) with C→N-masked `insert_left` matching, optional `require_c_all_t` conversion filter, demux mC→T QC stats, and conversion columns in [`scripts-v2/summary.py`](../../scripts-v2/summary.py). Configs: `workflow/dbit_taps_v2_test.json`, `workflow/dbit_taps_v2_slurm.json`. Unchanged stages reuse `scripts/*.py`.
- **mbias / call host subsample:** after fractional `samtools view -s` subsampling, host BAM preparation in [`scripts/host_subsample_bam.py`](../../scripts/host_subsample_bam.py) now caps alignment records at a fixed internal limit (`HOST_SUBSAMPLE_MAX_READS`, default 10M; not exposed in CLI or workflow JSON). TAPS and EMSeq `mbias` log `host_subsample_max_reads` for reproducibility.
- **TAPS CH calling:** TAPS `call` now supports `call_context_mode` values `cg`, `ch`, and `both`, with a dedicated CH caller at `scripts/methy_caller_CH.py` (later removed in 3.0.0). `cg` remains the default behavior.
- **TAPS CH outputs:** `.CH.cov` outputs now cover both reference-strand CH cytosines and opposite-strand CH cytosines represented on the reference strand as `DGN` anchors, instead of counting only reference-strand `CA/CC/CT` sites.
- **TAPS CH output contract:** `.CH.cov` rows now append `context` and `strand` columns. `context` is one of `CA`, `CC`, or `CT`, and `strand` is `+` or `-` to indicate which genomic strand contains the target cytosine.

## 2.0.0

- **make_cmd input checks:** TAPS (`scripts/make_cmd.py`), EMSeq (`scripts-emseq/make_cmd.py`), and RNA (`scripts-rna/make_cmd.py`) now verify that referenced files and directories exist before writing `.sh` / `.sbatch` scripts (FASTQs, indexes, references, whitelists, stage wrapper scripts, optional methscan paths, and—unless `--skip-workdir-input-checks` is set—expected prior-stage outputs under the sample work directory). For `--stage all`, TAPS/EMSeq run the same checks in the parent process first (with workdir intermediates skipped), then each per-stage subprocess; subprocess failures print stdout/stderr before exit. Shared helpers live in `scripts/workflow_input_checks.py`. See [config reference](../developers/config-reference.md).
- **BREAKING:** methscan outputs moved from `work/<sample>/coverage/` to `work/<sample>/methscan/` (`compact/`, `filter/`, `matrix/`, and default `TSS_profile.csv` / `VMRs.bed` at `methscan/` root). `methscan_profile_prepared_subdir` values are now `compact` or `filter` (not `host_prepare`). See [contracts](../developers/contracts.md#methscan-optional-stages-methscan_).
- added `scripts/methscan_run.py` with optional stages `methscan_filter`, `methscan_profile`, `methscan_smooth`, `methscan_scan`, `methscan_matrix`, and `methscan_all`
- fixed workflow JSON `methscan_smooth_use_weights` / `methscan_matrix_sparse` being ignored by TAPS and EMSeq `make_cmd` (CLI `store_true` defaults no longer override config when flags are omitted)
- remove `methscan_prepare` `chunksize` from workflow config

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
