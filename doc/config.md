# Config Guide

配置文件 `workflow/*.json` 说明。运行命令见 `doc/commands.md`，EMSeq 见 `doc/emseq.md`。

## 配置流程

1. 复制 `workflow/dbit_taps_test.json` 为新配置文件
2. 修改必填字段
3. 执行 `--dry-run` 验证路径解析
4. 按需调整线程数、并行度与 Slurm 资源

## 必填字段

- `sample_id`：样本名，决定 `work/<sample>/` 目录
- `r1`、`r2`：原始 FASTQ 路径
- `barcode1_whitelist`、`barcode2_whitelist`：barcode 白名单
- `bwa_index`：host 比对索引
- `call_reference_file`：calling 参考序列
- `spike_in_index`：可选，配置 spike-in 时需显式设置

Slurm 环境需额外配置各 stage 资源：

- `slurm.<stage>.partition`
- `slurm.<stage>.mem`
- `slurm.<stage>.cpus_per_task`

## 常用字段

### 运行控制

- `runner`：`local` 或 `slurm`
- `stage`：完整流程使用 `all`
- `work_root`：工作目录根路径
- `number_of_split_parts`：`fastp_split` chunk 数

### demux

- `linker_bc`：barcode2 与 barcode1 之间的 linker
- `insert_left`：insert 上游 anchor（等价于旧 `tn5` / `Tn5 mosaic end`）
- `linker_edit_distance`
- `barcode_hamming_distance`
- `gzip_level`

### align / pool / split

- `bwa_threads`
- `samtools_threads`
- `host_sort_mem`
- `split_barcodes`
- `split_smoke`：`true` 时最多输出 16 个非空 spot

### mbias / call / saturation / summary / aggregate（实验）

- `mbias_mode`：`spike`、`host`、`all`
- `mbias_host_subsample_fraction`
- `call_mode`：`spike`、`host`、`all`
- `call_r1_left_trimming` / `call_r1_right_trimming`
- `call_r2_left_trimming` / `call_r2_right_trimming`
- `saturation_reads_threshold`
- `aggregate_script`：可选，覆盖默认 `scripts/aggregate.py`
- `aggregate_sort_mem`：可选，传给 `scripts/aggregate.py` 中 GNU `sort -S` 的内存上限（默认 `8G`）

`call_r1_*_trimming` 与 `call_r2_*_trimming` 按原始 read 两端 cycle 定义，reverse-strand 比对保持 left/right 端语义。

## spike_in_index 格式

支持两种格式：

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

## EMSeq 配置差异

EMSeq 从 `workflow/dbit_emseq_test.json` 复制配置，主要差异：

- `fastp_split`：`sample_id`、`r1`、`r2`、`work_root`、`fastp_threads`、`number_of_split_parts`、`fastp_bin`
- `demux_extract_bc`：`barcode1_whitelist`、`barcode2_whitelist`、`linker_bc`、`insert_left`
- `align`：`biscuit_reference`、`biscuit_threads`、`biscuit_batch_size`、`biscuit_bin`、`sinto_bin`、`samtools_bin`
- **`split` 必填**：`split_barcodes`
- **`call` 必填**：`call_reference_file`、`call_jobs`
- **`call` 可选**：`call_left_trimming` / `call_right_trimming` → `biscuit pileup -5` / `-3`（每条 read 的 5'/3' 端最小距离；与 TAPS 的 `call_r1_*` / `call_r2_*` 不同，**不能**按 R1/R2 分别设置）。省略时不传对应 flag，使用 biscuit 内置默认。
- demux 可调参数：`linker_edit_distance`、`barcode_hamming_distance`、`gzip_level`
- Slurm 按 stage 配置：`slurm.split.split_bams`、`slurm.mbias.host`、`slurm.call.host`、`slurm.saturation`、`slurm.summary`；实验性 `aggregate`：`slurm.aggregate`、`aggregate_script`

可选字段详见 `doc/emseq.md`。

## 配置验证

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner local \
  --dry-run
```

输出 stage 顺序为 `fastp_split -> ... -> summary` 且无路径错误，表示配置有效。
