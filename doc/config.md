# Config Guide

本页说明 `workflow/*.json` 的常用配置，重点回答“首次跑自己的数据时哪些字段必须先改”。运行命令请看 `doc/commands.md`，EMSeq 独立入口请看 `doc/emseq.md`。

## 推荐做法

1. 复制 `workflow/dbit_taps_test.json` 到你自己的配置文件。
2. 先只修改“最小必改字段”。
3. 先执行 `--dry-run`，确认命令展开和路径解析正确。
4. 跑通后再调线程数、并行度和 Slurm 资源。

## 最小必改字段

下面字段通常必须先改：

- `sample_id`：样本名，会决定 `work/<sample>/` 目录名
- `r1`、`r2`：原始 FASTQ 路径
- `barcode1_whitelist`、`barcode2_whitelist`：barcode 白名单
- `bwa_index`：host 比对索引
- `call_reference_file`：calling 用参考
- `spike_in_index`：可选；有 spike-in 时建议显式设置

如果你跑 `slurm`，还需要确认每个 stage 的资源字段：

- `slurm.<stage>.partition`
- `slurm.<stage>.mem`
- `slurm.<stage>.cpus_per_task`

## 常改字段

### 运行控制

- `runner`：`local` 或 `slurm`
- `stage`：建议完整流程显式传 `all`
- `work_root`：工作目录根路径
- `number_of_split_parts`：`fastp_split` 的 chunk 数

### demux

- `linker1`、`linker2`、`tn5`
- `linker_edit_distance`
- `barcode_hamming_distance`
- `gzip_level`

### align / pool / split

- `bwa_threads`
- `samtools_threads`
- `host_sort_mem`
- `split_barcodes`
- `split_smoke`：`true` 时最多输出 16 个非空 spot

### mbias / call / saturation / summary

- `mbias_mode`：`spike`、`host`、`all`
- `mbias_host_subsample_fraction`
- `call_mode`：`spike`、`host`、`all`
- `call_r1_left_trimming` / `call_r1_right_trimming`
- `call_r2_left_trimming` / `call_r2_right_trimming`
- `saturation_reads_threshold`

`call_r1_*_trimming` 和 `call_r2_*_trimming` 按原始 read 两端的 cycle 定义生效，对 reverse-strand 比对同样保持 left/right 端语义。

## spike_in_index 写法

支持两种写法：

1. JSON 对象：

```json
"spike_in_index": {
  "lambda": "/path/to/lambda.fa",
  "puc19": "/path/to/puc19.fa"
}
```

2. `NAME=INDEX` 列表：

```text
--spike-in-index lambda=/path/to/lambda.fa
--spike-in-index puc19=/path/to/puc19.fa
```

## EMSeq 入口的配置差异

EMSeq 当前支持 `fastp_split -> demux_extract_bc -> align`，应从 `workflow/dbit_emseq_test.json` 复制配置，并重点关注：

- `fastp_split` 最小字段：`sample_id`、`r1`、`r2`、`work_root`、`fastp_threads`、`number_of_split_parts`、`fastp_bin`
- `demux_extract_bc` 额外字段：`barcode1_whitelist`、`barcode2_whitelist`、`linker1`、`linker2`、`tn5`
- `align` 额外字段：`biscuit_reference`、`biscuit_threads`、`biscuit_batch_size`、`biscuit_bin`、`sinto_bin`、`samtools_bin`
- demux 的可调参数与 `scripts-emseq/extract_bc.py` 对齐：`linker_edit_distance`、`barcode_hamming_distance`、`gzip_level`
- Slurm 场景按 stage 分别配置 `slurm.fastp_split.*`、`slurm.demux_extract_bc.*` 与 `slurm.align.*`

`workflow/dbit_emseq_test.json` 中若出现 `pool`、`split`、`mbias`、`call`、`saturation`、`summary` 等后续阶段字段，应视为预留字段，不表示当前已经支持这些 stage。

## 最小检查命令

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner local \
  --dry-run
```

如果你看到 stage 顺序为
`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> saturation -> summary`，且无路径错误，说明配置基本可用。
