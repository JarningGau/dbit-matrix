# EMSeq Guide

本页是 `EMSeq` 独立入口的用户说明，覆盖当前 MVP 范围。维护者回归检查请看 `TEST-emseq.md`。

## 当前支持范围

- 当前状态：试验性 / MVP
- 入口脚本：`scripts-emseq/make_cmd.py`
- 当前支持 stage：`fastp_split -> demux_extract_bc -> align`
- 当前不包含：`pool -> split -> mbias -> call -> saturation -> summary`

## 与 TAPS 的主要差异

- `EMSeq` 使用独立入口，不通过 `scripts/make_cmd.py`
- `demux_extract_bc` 使用 `scripts-emseq/extract_bc.py`
- `align` 使用 `scripts-emseq/aligner.py`
- EMSeq 假定的 `R1` 结构为：
  `linker1-barcodeB-linker2-barcodeA-others(15 bp)-Tn5-insert`
- `align` 使用 `biscuit align` 进行比对，并将 host 结果整理为带 `CB` tag 的 BAM
- demux 和 align 的外部参数风格与 TAPS 对齐，当前仍通过 `workflow/*.json` 配置

## 最小配置

建议从 `workflow/dbit_emseq_test.json` 复制，并至少确认以下字段：

- `sample_id`
- `r1`
- `r2`
- `work_root`
- `fastp_threads`
- `number_of_split_parts`
- `fastp_bin`
- `barcode1_whitelist`
- `barcode2_whitelist`
- `linker1`
- `linker2`
- `tn5`
- `biscuit_reference`
- `biscuit_threads`
- `biscuit_batch_size`
- `biscuit_bin`
- `sinto_bin`
- `samtools_bin`

可调的 demux 参数：

- `linker_edit_distance`
- `barcode_hamming_distance`
- `gzip_level`

可调的 align 参数：

- `spike_in_index`

`workflow/dbit_emseq_test.json` 中若出现 `pool`、`split`、`mbias`、`call`、`saturation`、`summary` 等后续阶段字段，应视为预留字段，不表示当前已经支持这些 stage。

## 首次运行命令

先检查 `fastp_split`：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage fastp_split \
  --runner local \
  --dry-run
```

再检查 `demux_extract_bc`：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage demux_extract_bc \
  --runner local \
  --dry-run
```

如需提交到 Slurm，把 `--runner local` 改为 `--runner slurm`。

最后检查 `align`：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage align \
  --runner local \
  --dry-run
```

## 预期输出

`fastp_split` 后：

- `work/<sample>/shard_fastq/*.R1.fq.gz`
- `work/<sample>/shard_fastq/*.R2.fq.gz`
- `work/<sample>/shard_fastq/fastp.html`
- `work/<sample>/shard_fastq/fastp.json`

`demux_extract_bc` 后：

- `work/<sample>/demux/*.R1.demux.fq.gz`
- `work/<sample>/demux/*.R2.demux.fq.gz`
- `work/<sample>/demux/*.R1.spike-in.fq.gz`
- `work/<sample>/demux/*.R2.spike-in.fq.gz`
- `work/<sample>/demux/*.stats.json`

`align` 后：

- `work/<sample>/align_shards/*.cb.bam`
- `work/<sample>/align_shards/*.<spike_name>.bam`
