# Stage Contracts

This document is the canonical stage input/output reference for developers.

## TAPS And EMSeq Mainline

Main stage order:

`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> saturation -> summary`

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

- host: `pooled/pooled.byCB.bam`
- spike-in: `pooled/pooled.<spike_name>.sorted.bam`
- spike-in index: `pooled/pooled.<spike_name>.sorted.bam.bai`

### `split`

Purpose: split the pooled host BAM into per-spot BAMs.

Input:

- `pooled/pooled.byCB.bam`

Outputs:

- `split_bams/<X_index>/<X_index>_<Y_index>.bam`
- `split_bams/per_spot_read_counts.tsv`

Contract:

- spot parsing uses `CB:Z:<x>+<y>`
- smoke mode writes at most 16 non-empty spot BAMs
- in Slurm mode, `split_bams` and `sort` are separate jobs

### `split` sort substep

Purpose: sort and index per-spot BAMs.

Outputs:

- `split_bams/**/*.sorted.bam`
- `split_bams/**/*.sorted.bam.bai`

Contract:

- this is a substep of `split`, not a top-level stage

### `mbias`

Purpose: write M-bias QC outputs and prepare a host subsampled BAM when needed.

Inputs:

- host: `pooled/pooled.byCB.bam`
- spike-in: `pooled/pooled.<spike_name>.sorted.bam`

Outputs:

- `qc/mbias/host.subsampled.sorted.bam`
- `qc/mbias/host.subsampled.sorted.bam.bai`
- `qc/mbias/host.mbias.tsv`
- `qc/mbias/host.mbias.png`
- `qc/mbias/<spike_name>.mbias.tsv`
- `qc/mbias/<spike_name>.mbias.png`

Workflow-specific note:

- EMSeq uses `scripts-emseq/mbias.py`
- TAPS uses `scripts/mbias.py`

### `call`

Purpose: generate host per-spot, host mitochondrial aggregate, and spike-in methylation outputs.

Inputs:

- host per-spot: `split_bams/**/*.sorted.bam`
- host mito fallback source: `qc/mbias/host.subsampled.sorted.bam`
- spike-in: `pooled/pooled.<spike_name>.sorted.bam`

Outputs:

- `coverage/host/<X_index>/<X_index>_<Y_index>.CG.cov`
- `coverage/host/<X_index>/<X_index>_<Y_index>.CH.cov` (optional)
- `coverage/host_mito.CG.cov`
- `coverage/host_mito.CH.cov` (optional)
- `coverage/<spike_name>.CG.cov`
- `coverage/<spike_name>.CH.cov` (optional)

Contract:

- host outputs are per-spot
- `host_mito.CG.cov` is a single aggregate output
- `host_mito.CH.cov` is a single aggregate output when CH calling is enabled
- `.CG.cov` outputs contain covered CpG sites only
- `.CH.cov` outputs are headerless, include covered CH sites only, and add trailing `context` and `strand` columns
- `context` is one of `CA`, `CC`, or `CT`
- `strand` is `+` or `-` and indicates which genomic strand contains the target cytosine
- CH outputs are controlled by `call_context_mode`

### `saturation`

Purpose: estimate CpG saturation from host coverage and per-spot read counts.

Inputs:

- `coverage/host/**/*.CG.cov`
- `split_bams/per_spot_read_counts.tsv`

Outputs:

- `qc/saturation/saturation_curve.png`
- `qc/saturation/saturation_summary.tsv`

### `summary`

Purpose: generate per-spot and per-sample reports.

Inputs include:

- host coverage
- host mitochondrial coverage
- spike-in coverage
- split read counts
- fastp statistics
- demux statistics
- pooled BAMs
- saturation summary

Outputs:

- `summary/per_spot_summary.tsv`
- `summary/sample_summary.tsv`
- `summary/reads_heatmap.png`
- `summary/cpg_site_count_heatmap.png`
- `summary/mean_methylation_heatmap.png`

Contract:

- `per_spot_summary.tsv` is one row per spot
- `sample_summary.tsv` is one row per sample
- missing upstream inputs should preserve fixed columns and write `NA`
- `per_spot_summary.tsv` keeps fixed columns for `X_index`, `Y_index`, `spot`, `mean_methylation`, `cpg_site_count`, and `reads`
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
- `mbias` and `call` have EMSeq-specific chemistry and trimming behavior, but the output contracts above remain stable

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
