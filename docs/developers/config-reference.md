# Config Reference

This document defines the workflow configuration surface for developers.

## General Rules

- workflow parameters should be stored in `workflow/*.json`
- all stages must support `--dry-run`
- stage names in docs should match the CLI values exactly
- Slurm resources should be configured per stage

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
- `methscan_prepare_script`
- `methscan_pixi_manifest`

These stages are optional and are not part of `--stage all`.

## Slurm Structure

Slurm resources should remain stage-specific. Examples:

- `slurm.mbias.host`
- `slurm.mbias.spike`
- `slurm.call.host`
- `slurm.call.spike`
- `slurm.saturation`
- `slurm.summary`

Do not collapse unrelated stages into one shared resource block.

For RNA, stage logs default to `work/<sample>/logs/` and can be overridden through workflow or CLI log settings.
