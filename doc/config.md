# Config Guide

本页说明 `workflow/*.json` 的常用配置。目标是让首次用户知道哪些字段必须先改，哪些字段可以后续再调。

## 推荐做法

1. 复制 `workflow/dbit_taps_test.json` 到你自己的配置文件。
2. 先只修改“最小必改字段”。
3. 先执行 `--dry-run`，确认命令展开和路径解析正确。
4. 跑通后再调性能参数（线程、并行度、slurm 资源）。

## 最小必改字段

下面字段是首次跑自己的数据时通常必须修改的：

- `sample_id`：样本名，会决定 `work/<sample>/` 目录名
- `r1`、`r2`：原始 FASTQ 路径
- `barcode1_whitelist`、`barcode2_whitelist`：barcode 白名单
- `bwa_index`：host 比对索引
- `call_reference_file`：calling 用参考
- `spike_in_index`：可选；有 spike-in 时建议显式设置

如果你跑 Slurm，还需要确认：

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
- `split_smoke`（true 时最多输出 16 个非空 spot）

### mbias / call / summary

- `mbias_mode`：`spike`、`host`、`all`
- `mbias_host_subsample_fraction`
- `call_mode`：`spike`、`host`、`all`
- `call_r1_left_trimming` / `call_r1_right_trimming`
- `call_r2_left_trimming` / `call_r2_right_trimming`
  以上 trimming 参数按原始 read 两端的 cycle 定义生效，对 reverse-strand 比对同样保持 left/right 端语义。

## spike_in_index 写法

支持两种写法：

1) JSON 对象（推荐，最清晰）

```json
"spike_in_index": {
  "lambda": "/path/to/lambda.fa",
  "puc19": "/path/to/puc19.fa"
}
```

2) `NAME=INDEX` 列表（CLI 传参场景）

```text
--spike-in-index lambda=/path/to/lambda.fa
--spike-in-index puc19=/path/to/puc19.fa
```

## 最小检查命令

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner local \
  --dry-run
```

如果你看到 stage 顺序为  
`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> summary`，且无路径错误，说明配置基本可用。
