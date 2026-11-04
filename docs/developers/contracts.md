# Stage Contracts

This document is the canonical stage input/output reference for developers.

## TAPS And EMSeq Mainline

Main stage order:

`fastp_split -> demux_extract_bc -> align -> pool -> mbias -> call -> saturation -> summary`

### `fastp_split`

Purpose: quality control and chunking of paired FASTQ input.

Inputs:

- raw `R1 FASTQ`
- raw `R2 FASTQ`

Outputs:

- `work/<sample>/shard_fastq/*.R1.fq.gz`
- `work/<sample>/shard_fastq/*.R2.fq.gz`
- `work/<sample>/shard_fastq/fastp.html`
- `work/<sample>/shard_fastq/fastp.json`

### `demux_extract_bc`

Purpose: extract DBiT barcodes, separate matched reads and spike-in reads, and write demux statistics.

Inputs:

- `shard_fastq/*.R1.fq.gz`
- `shard_fastq/*.R2.fq.gz`

Outputs:

- `demux/<chunk>.R1.demux.fq.gz`
- `demux/<chunk>.R2.demux.fq.gz`
- `demux/<chunk>.R1.spike-in.fq.gz`
- `demux/<chunk>.R2.spike-in.fq.gz`
- `demux/<chunk>.stats.json`

Contract:

- stats must include retention information and reject reasons
- matched and spike-in reads must be separated
- FASTQ read names carry CB as concatenated `<bcB><bcA>` (no `+`): `@<bcB><bcA>:<original>`
- barcode1 = barcodeA (X); barcode2 = barcodeB (Y); physical R1 order remains bcB–linker–bcA

#### TAPS v2 demux (`scripts-v2/extract_bc.py`)

TAPS v2 uses the same output paths and mainline stage order. The demux implementation additionally:

- locates `insert_left` (the methylated-linker on R1) with C→N masking so converted anchors still match
- always records mC→T conversion QC in `stats.json` (`spike_mCtoT_*`, `spike_scored_reads`, `scored_fully_converted`, …)
- optionally rejects reads that fail an all-T conversion filter (`--require-c-all-t` / `require_c_all_t`; empty = off)

Terminology: the methylated-linker is the 5mC-bearing `insert_left` conversion QC motif. It is not a reference spike-in (lambda/pUC19). Rejected reads still use the existing `*.spike-in.fq.gz` contract names.

### `align`

Purpose: align demultiplexed reads to spike-in and host references.

Inputs:

- host: `demux/*.demux.fq.gz`
- spike-in: `demux/*.spike-in.fq.gz`

Outputs:

- host: `align_shards/<chunk>.cb.bam`
- spike-in: `align_shards/<chunk>.<spike_name>.bam`

Contract:

- host input must come from `*.demux.fq.gz`
- spike-in input must come from `*.spike-in.fq.gz`
- the stage runs spike-in alignment before host alignment
- `spike_in_index` must support either a JSON object or a `NAME=INDEX` list

### `pool`

Purpose: merge host and spike-in BAMs across chunks.

Inputs:

- host: `align_shards/*.cb.bam`
- spike-in: `align_shards/*.<spike_name>.bam`

Outputs:

- host: `pooled/pooled.host.sorted.bam`
- host index: `pooled/pooled.host.sorted.bam.bai` (or `.csi`)
- spike-in: `pooled/pooled.<spike_name>.sorted.bam`
- spike-in index: `pooled/pooled.<spike_name>.sorted.bam.bai`

Contract:

- host is coordinate-sorted (no `-t CB`) and indexed
- spike-in `pooled/pooled.<spike_name>.sorted.bam` contains mapped reads only (`samtools view -F 4` before sort; secondary/supplementary mapped reads are kept)

### `split` (legacy script, not a workflow stage)

`scripts/split_bams.py` remains in the repository. It is not in `--stage all` and is not a `STAGE_CHOICES` value. When used, it writes `split_bams/<cb[0:3]>/<cb>.bam` and `split_bams/per_spot_read_counts.tsv` (`barcode`, `reads`). Per-spot BAM trees are not a mainline output.

### `mbias`

Purpose: write M-bias QC outputs and prepare a host subsampled BAM when needed.

Inputs:

- host: `pooled/pooled.host.sorted.bam`
- spike-in: `pooled/pooled.<spike_name>.sorted.bam`

Outputs:

- `qc/mbias/host.subsampled.sorted.bam`
- `qc/mbias/host.subsampled.sorted.bam.bai`
- `qc/mbias/host.mbias.tsv`
- `qc/mbias/host.mbias.png`
- `qc/mbias/<spike_name>.mbias.tsv`
- `qc/mbias/<spike_name>.mbias.png`

Host counting note:

- Host M-bias CG counts use the nuclear chromosome whitelist from `call_chromosomes`
- Mitochondrial contigs such as `chrM` are excluded from `host.mbias.tsv`
- Only cycles with `coverage > 500` are included in TSV output and plots

Workflow-specific note:

- EMSeq uses `scripts-emseq/mbias.py`
- TAPS uses `scripts/mbias.py`

### `call`

Purpose: generate host per-spot, host mitochondrial aggregate, and spike-in methylation outputs.

Inputs:

- host: `pooled/pooled.host.sorted.bam` (indexed) plus `split_barcodes`
- host mito fallback source: `qc/mbias/host.subsampled.sorted.bam`
- spike-in: `pooled/pooled.<spike_name>.sorted.bam`

Outputs:

- `coverage/per_spot_read_counts.tsv` (`barcode`, `reads`; `barcode` = CB `<bcB><bcA>`)
- `coverage/host/<cb[0:3]>/<cb>.CG.cov`
- `coverage/host/<cb[0:3]>/<cb>.CH.cov` (optional)
- `coverage/host_mito.CG.cov`
- `coverage/host_mito.CH.cov` (optional)
- `coverage/<spike_name>.CG.cov`
- `coverage/<spike_name>.CH.cov` (optional)

Contract:

- `split_barcodes` is the CB whitelist for `call`, not a split-stage input
- host calling does not read per-spot BAMs
- host outputs are per-spot; spot identity is the BAM `CB` tag (`<bcB><bcA>`)
- host cov paths use the first three bases of CB as the subdirectory (64 possible prefix dirs)
- `call` writes `coverage/per_spot_read_counts.tsv` from a full BAM scan of whitelist-assigned reads
- smoke is `call_smoke`: write `.cov` only for the 16 spots with highest `reads` (ties broken by CB ascending)
- `host_mito.CG.cov` is a single aggregate output
- `host_mito.CH.cov` is a single aggregate output when CH calling is enabled
- `.CG.cov` outputs contain covered CpG sites only
- `.CH.cov` outputs are headerless, include covered CH sites only, and add trailing `context` and `strand` columns
- `context` is one of `CA`, `CC`, or `CT`
- `strand` is `+` or `-` and indicates which genomic strand contains the target cytosine
- CH outputs are controlled by `call_context_mode`
- EMSeq host/mito/spike calling uses pooled pysam (`scripts-emseq/methy_caller.py` / `scripts-emseq/call.py`), not biscuit pileup or VCF intermediates
- EMSeq `.CG.cov` positions are pysam 0-based; `pct` has two decimal places
- EMSeq trimming uses `call_trimming` like TAPS (see [config reference](config-reference.md#emseq-fields))

### `saturation`

Purpose: estimate CpG saturation from host coverage and per-spot read counts.

Inputs:

- `coverage/host/**/*.CG.cov`
- `coverage/per_spot_read_counts.tsv`
- `shard_fastq/fastp.json` (optional; supplies total sequencing depth for the Gbp x-axis. If absent, the plot falls back to a raw coverage-fraction x-axis.)

HQ spot filter: `reads > reads_threshold` and a matching `*.CG.cov`. `reads_threshold` comes from `--reads-threshold` / config `saturation_reads_threshold`: `"auto"` (default) or a positive number. Auto sorts barcodes by reads descending, drops `reads <= 10`, then takes the reads value at the steepest slope on `(rank, log10(reads+1))` (smoothed window 51). See [config reference](config-reference.md#taps-fields).

Outputs:

- `qc/saturation/saturation_curve.png`
- `qc/saturation/reads_knee.png` (barcode rank vs `log10(reads+1)` with the resolved threshold marked)
- `qc/saturation/saturation_summary.tsv` (columns: `sample_id`, `observed_median_unique_cpgs`, `theoretical_max_median_unique_cpgs`, `predicted_median_unique_cpgs_at_2x`, `saturation_rate`, `extrapolation_model`, `hq_spot_count`, `reads_threshold`, `reads_threshold_mode`)

Extrapolation model: fits both a through-origin linear model and an exponential saturation curve. When the linear fit reaches `--linear-r2-threshold` (default `0.99`) the curve is treated as undersaturated, `extrapolation_model=linear` is reported, and `saturation_rate`/theoretical max are `NA`; otherwise the saturation curve is used.

### `summary`

Purpose: generate per-spot and per-sample reports.

Inputs include:

- host coverage
- host mitochondrial coverage
- spike-in coverage
- per-spot read counts (`coverage/per_spot_read_counts.tsv`, written by `call`)
- `barcode_positions` TSV (`barcode`, `array_row`, `array_col`) where `barcode` = `<bcB><bcA>`
- fastp statistics
- demux statistics
- pooled BAMs
- saturation summary

Outputs:

- `summary/per_spot_summary.tsv`
- `summary/sample_summary.tsv`
- `summary/reads_heatmap.png`
- `summary/cpg_site_count_heatmap.png`
- `summary/mean_5mCG_heatmap.png` (colorbar 50–80%)
- `summary/mean_5mCH_heatmap.png` (colorbar 0–10%)

Contract:

- `per_spot_summary.tsv` is one row per spot
- `sample_summary.tsv` is one row per sample
- missing upstream inputs should preserve fixed columns and write `NA`
- `per_spot_summary.tsv` fixed columns: `barcode`, `array_row`, `array_col`, `mean_5mCG`, `mean_5mCH`, `cpg_site_count`, `ch_site_count`, `reads`
- `mean_5mCG` / `mean_5mCH` are unweighted means of site `%` from per-spot `*.CG.cov` / `*.CH.cov`; missing or empty cov writes `NA` for the matching rate and site count
- summary `barcode` is `<bcB><bcA>` (same CB as cov stem / counts; no `host` suffix)
- `array_row` / `array_col` are copied from `barcode_positions` by joining on CB; missing positions write `NA`
- heatmaps plot metrics on `array_col` (x) and `array_row` (y), 0-based inclusive
- all summary heatmaps place `(array_col=0, array_row=0)` at the top-right (`array_col` increases leftward, `array_row` increases downward)
- `mean_5mCG_heatmap.png` / `mean_5mCH_heatmap.png` plot only HQ spots (`reads > reads_threshold` from `qc/saturation/saturation_summary.tsv`); non-HQ and missing values stay gray (`#f2f2f2`). Grid extent uses all spots with valid `array_row`/`array_col`. `reads_heatmap.png` and `cpg_site_count_heatmap.png` are unchanged. If `reads_threshold` is missing, methylation heatmaps fall back to all spots with a warning.
- `X_index`, `Y_index`, and `spot` (`XX_YY`) are removed
- `sample_summary.tsv` methylation fields are context-split: `host_spot_mean_5mCG` / `host_spot_mean_5mCH` (weighted by `cpg_site_count` / `ch_site_count`), `host_mito_mean_5mCG` / `host_mito_mean_5mCH`, and `{spike}_mean_5mCG` / `{spike}_mean_5mCH`; `host_spot_median_cpg_sites` is unchanged
- `sample_summary.tsv` includes `host_valid_reads`, plus the mito/nuclear split `host_valid_reads_mito` and `host_valid_reads_nuclear` (same valid-read filters; mito contigs from `--mito-chromosomes`, default `chrM`; nuclear is the complement). When BAM metrics succeed: `host_valid_reads_mito + host_valid_reads_nuclear == host_valid_reads`
- TAPS v2 (`scripts-v2/summary.py`) additionally surfaces methylated-linker conversion QC columns in `sample_summary.tsv`: `spike_scored_reads`, `scored_fully_converted`, `spike_mean_conversion_rate`, `spike_fraction_all_T`, `spike_conversion_filtered_rejects`

## Optional TAPS And EMSeq Stages

### `aggregate`

Purpose: flatten host per-spot `.CG.cov` rows into aggregated TSV outputs.

Input:

- `coverage/host/**/*.CG.cov`

Outputs:

- `coverage/aggregated_cg_by_id.tsv`
- `coverage/aggregated_cg_by_pos.tsv`

Contract:

- outputs are headerless TSV files
- both outputs keep the same six logical fields: `id`, `chr`, `start`, `end`, `mC`, `C`
- `aggregated_cg_by_id.tsv` is sorted by `id`, then position
- `aggregated_cg_by_pos.tsv` is sorted by position, then `id`

### Methscan optional stages (`methscan_*`)

All methscan steps are **host-only**, **explicit** (not part of `--stage all`), and run `scripts/methscan_run.py` in the `envs/methscan` pixi workspace unless `methscan_pixi_manifest` overrides it. Per-cell coverage **input** stays under `work/<sample>/coverage/host/`; all **methscan outputs** live under `work/<sample>/methscan/`.

#### `methscan_prepare`

Purpose: run `methscan prepare` on per-spot Bismark-like coverage.

Input:

- `coverage/host/**/*.CG.cov`

Output:

- `methscan/compact/`

#### `methscan_filter`

Purpose: run `methscan filter` from prepared data to filtered cells.

Input:

- `methscan/compact/`

Output:

- `methscan/filter/`

#### `methscan_profile`

Purpose: run `methscan profile` for mean methylation over regions (e.g. TSS).

Input:

- sorted regions `.bed` (`methscan_tss_bed` in workflow config or CLI)
- prepared directory: default `methscan/compact/`; optional `methscan/filter/` via `methscan_profile_prepared_subdir` (`compact` or `filter`)

Output:

- default `methscan/TSS_profile.csv` (override with `methscan_profile_csv`)

#### `methscan_smooth`

Purpose: run `methscan smooth` on filtered data.

Input:

- `methscan/filter/`

Output:

- smoothed values under `methscan/filter/smoothed/` (methscan layout)

#### `methscan_scan`

Purpose: run `methscan scan` for VMRs.

Input:

- `methscan/filter/` (includes smoothed outputs from `methscan smooth`)

Output:

- default `methscan/VMRs.bed` (override with `methscan_vmrs_bed`)

#### `methscan_matrix`

Purpose: run `methscan matrix` for region × cell tables.

Input:

- VMRs `.bed` (default `methscan/VMRs.bed`)
- `methscan/filter/`

Output:

- default output directory `methscan/matrix/` (override with `methscan_matrix_prefix`)

#### `methscan_all`

Purpose: run `prepare`, `filter`, `profile`, `smooth`, `scan`, and `matrix` **in order** in one generated script (same effect as running those stages sequentially). Requires `methscan_tss_bed` for the profile step.

## EMSeq-Specific Additions

- `--stage all` covers only the fixed mainline
- `aggregate` and the methscan optional stages are not part of `all`
- `mbias` has EMSeq-specific chemistry (`scripts-emseq/mbias.py`)
- `call` polarity is inverted relative to TAPS (retained base = methylated); `.cov` paths still follow the contract above
- EMSeq uses `call_trimming` like TAPS; biscuit remains an align-stage tool, not a call-stage tool

## RNA Current Contracts

### `demux_extract_bc`

Inputs:

- raw `r1`
- raw `r2`

Outputs:

- `demux/<sample>.R1.clean.fq.gz`
- `demux/<sample>.R2.clean.fq.gz`
- `demux/<sample>.stats.json`

Contract:

- RNA does not use chunk inputs
- in Slurm mode, default stage logs are written under `work/<sample>/logs/`

### `align`

Inputs:

- `demux/<sample>.R1.clean.fq.gz`
- `demux/<sample>.R2.clean.fq.gz`

Outputs:

- `solo/`
- `solo/star.Solo.out/Gene/raw/barcodes_pos.tsv`
- `solo/star.Solo.out/Gene/filtered/barcodes_pos.tsv`

Contract:

- RNA currently produces STARsolo matrix outputs
- RNA does not produce methylation BAM or `.CG.cov` outputs
- in Slurm mode, default stage logs are written under `work/<sample>/logs/`
