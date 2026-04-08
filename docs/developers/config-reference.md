# Config Reference

This document defines the workflow configuration surface for developers.

## General Rules

- workflow parameters should be stored in `workflow/*.json`
- all stages must support `--dry-run`
- stage names in docs should match the CLI values exactly
- Slurm resources should be configured per stage
- TAPS (`scripts/make_cmd.py`), EMSeq (`scripts-emseq/make_cmd.py`), and RNA (`scripts-rna/make_cmd.py`) validate that referenced input paths exist before emitting `.sh` / `.sbatch` scripts. Use `--skip-workdir-input-checks` when generating a full pipeline driver (`--stage all`) or a downstream stage before upstream outputs exist under the sample work directory; `--stage all` passes this flag to each per-stage subprocess automatically (TAPS/EMSeq). For `--stage all`, TAPS/EMSeq also run the same path checks in the parent process first (with workdir intermediates skipped), so invalid `r1`/`r2` or missing references fail immediately with a clear error before any subprocess runs.

## Common Fields

### Run Control

- `sample_id`: sample name used to build `work/<sample>/`
- `r1`, `r2`: input FASTQ paths
- `work_root`: working directory root
- `runner`: `local` or `slurm`
- `stage`: named stage or `all`
- `number_of_split_parts`: chunk count for `fastp_split` when the workflow uses chunking

### Demux

- `barcode1_whitelist`
- `barcode2_whitelist`
- `linker_bc`
- `insert_left`
- `linker_edit_distance`
- `barcode_hamming_distance`
- `gzip_level`

## TAPS Fields

Core fields:

- `bwa_index`
- `call_reference_file`
- `number_of_split_parts`

Common performance fields:

- `bwa_threads`
- `samtools_threads`
- `host_sort_mem`
- `split_barcodes`
- `split_smoke`

Calling fields:

- `call_mode`
- `call_r1_left_trimming`
- `call_r1_right_trimming`
- `call_r2_left_trimming`
- `call_r2_right_trimming`

Notes:

- TAPS trimming is defined separately for R1 and R2 ends

## EMSeq Fields

Core fields:

- `biscuit_reference`
- `split_barcodes`
- `call_reference_file`
- `call_jobs`

Common performance fields:

- `fastp_threads`
- `number_of_split_parts`
- `biscuit_threads`
- `biscuit_batch_size`
- `samtools_threads`
- `host_sort_mem`

Calling fields:

- `call_mode`
- `call_left_trimming`
- `call_right_trimming`
- `call_host_threads`
- `call_spike_threads`
- `call_host_subsample_fraction`
- `call_host_subsample_seed`

QC fields:

- `mbias_mode`
- `mbias_host_subsample_fraction`
- `mbias_max_cycle`
- `mbias_min_mapping_quality`

Notes:

- EMSeq trimming uses `call_left_trimming` and `call_right_trimming`
- EMSeq does not split trimming by R1 and R2

## RNA Fields

Core fields:

- `barcode1_whitelist`
- `barcode2_whitelist`
- `linker_bc`
- `umi_left`
- `umi_len`
- `star_genome_dir`
- `gtf`

## Spike-In Configuration

`spike_in_index` supports:

- a JSON object
- a `NAME=INDEX` list

Example object:

```json
"spike_in_index": {
  "lambda": "/path/to/lambda.fa",
  "puc19": "/path/to/puc19.fa"
}
```

Example repeated CLI form:

```bash
--spike-in-index lambda=/path/to/lambda.fa \
--spike-in-index puc19=/path/to/puc19.fa
```

## Optional Advanced Stages

TAPS and EMSeq may also use:

- `aggregate_script`
- `aggregate_sort_mem`
- `methscan_run_script` (default: `scripts/methscan_run.py`)
- `methscan_pixi_manifest`
- `methscan_prepare_chunksize` (default: `10000000`)
- `methscan_filter_min_sites` (default: `50000`)
- `methscan_tss_bed` (required for `methscan_profile` and `methscan_all`)
- `methscan_profile_strand_column` (default: `6`)
- `methscan_profile_prepared_subdir`: `compact` or `filter` (default: `compact`)
- `methscan_profile_csv` (optional path)
- `methscan_smooth_bandwidth` (optional)
- `methscan_smooth_use_weights` (boolean)
- `methscan_scan_threads` (default: `10`)
- `methscan_vmrs_bed` (optional path)
- `methscan_matrix_threads` (default: `10`)
- `methscan_matrix_sparse` (boolean)
- `methscan_matrix_prefix` (optional output directory)

For `methscan_smooth_use_weights` and `methscan_matrix_sparse`, omit the matching CLI flags so values from workflow JSON apply; passing `--methscan-smooth-use-weights` or `--methscan-matrix-sparse` on the command line still forces those options on.

These stages are optional and are not part of `--stage all`.

## Slurm Structure

Slurm resources should remain stage-specific. Examples:

- `slurm.mbias.host`
- `slurm.mbias.spike`
- `slurm.call.host`
- `slurm.call.spike`
- `slurm.saturation`
- `slurm.summary`
- `slurm.methscan_prepare`, `slurm.methscan_filter`, `slurm.methscan_profile`, `slurm.methscan_smooth`, `slurm.methscan_scan`, `slurm.methscan_matrix`, `slurm.methscan_all` (methscan optional stages)

Do not collapse unrelated stages into one shared resource block.

For RNA, stage logs default to `work/<sample>/logs/` and can be overridden through workflow or CLI log settings.
