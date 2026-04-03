# Outputs Guide

结果解读与排错指南。完整路径与 stage 契约见 `doc/stages.md`。TAPS 与 EMSeq 的 `saturation` 与 `summary` 产物路径一致。

## 检查顺序

1. `summary/sample_summary.tsv`：样本级总览
2. `summary/per_spot_summary.tsv`：spot 级甲基化与 reads
3. `summary/*.heatmap.png`：可视化
4. `qc/saturation/saturation_summary.tsv` 与 `qc/saturation/saturation_curve.png`：饱和度

## Summary 文件

### `per_spot_summary.tsv`

每 spot 一行，固定列：

- `X_index`、`Y_index`、`spot`
- `mean_methylation`、`cpg_site_count`、`reads`

### `sample_summary.tsv`

每样本一行；缺失输入保持固定列，对应位置写 `NA`。常见列：

- `raw_reads`、`barcoded_reads`、`barcoded_reads_rate`
- `host_mapped_reads`、`host_valid_reads`、`valid_reads_rate`
- `<spike_name>_mapped_reads`（配置 spike 时）
- `saturation_rate`：来自 `qc/saturation/saturation_summary.tsv`；未执行 `saturation` 或无结果时为 `NA`

reads 数值可能为千分位格式。

## Calling 与 QC 文件

### Host per-spot coverage

- `work/<sample>/coverage/host/<X_index>/<X_index>_<Y_index>.CG.cov`
- **EMSeq**：默认剔除线粒体 contig（默认 `chrM`，由 `call_mito_chromosomes` 配置）

### Host mito 聚合

- `work/<sample>/coverage/host_mito.CG.cov`
- **EMSeq**：优先由 `qc/mbias/host.subsampled.sorted.bam` 经 pileup 生成；若该 BAM 不存在，则从各 spot `coverage/host/**/*.CG.cov` 汇总线粒体位点（默认 `chrM`）
- 上游依赖：`call`；与 `mbias` 是否产出 host subsampled BAM 相关，见 `doc/stages.md`

### Spike 聚合（可选）

- `work/<sample>/coverage/<spike_name>.CG.cov`（配置 `spike_in_index` 且执行对应 calling）

### `aggregate` 扁平表（实验性，可选）

- `work/<sample>/coverage/aggregated_cg_by_id.tsv`、`aggregated_cg_by_pos.tsv`：由 TAPS 或 EMSeq 编排的 `aggregate` stage（`scripts/make_cmd.py` / `scripts-emseq/make_cmd.py`）从 `coverage/host/**/*.CG.cov` 生成；**无表头**，列顺序与含义见 `doc/stages.md` §11

### `methscan_prepare`（实验性，可选）

- `work/<sample>/coverage/host_prepare/`：由 `methscan prepare`（经 `scripts/methscan_run.py prepare`、pixi 工作区 `envs/methscan`）从 `coverage/host/**/*.CG.cov` 生成；仅 host，契约见 `doc/stages.md` §12

### Saturation

- `work/<sample>/qc/saturation/saturation_summary.tsv`
- `work/<sample>/qc/saturation/saturation_curve.png`
- `summary` 读取 `saturation_rate` 写入 `sample_summary.tsv`

### mbias QC

- `work/<sample>/qc/mbias/*.mbias.tsv`、`.mbias.png`
- 用于 M-bias 诊断；与 `call` 中 `host_mito` 是否优先用 subsampled BAM 相关，见 `doc/stages.md`

## 常见问题排查

| 现象 | 检查项 |
|------|--------|
| `reads` 极低或大面积为 0 | `demux` 保留率、`demux/*.stats.json`、`split` 输入 |
| spike-in 相关缺失 | `spike_in_index`、`align`/`pool`/`call_mode` |
| `host_mito.CG.cov` 缺失或异常 | `call` 是否成功；`qc/mbias/host.subsampled.sorted.bam` 是否存在 |
| `saturation_rate` 为 `NA` | 是否执行 `saturation`；`coverage/host/**/*.CG.cov` 与 `split_bams/per_spot_read_counts.tsv` 是否存在 |

## 相关文档

- `doc/stages.md`：各 stage 输入输出契约
- `doc/commands.md`：TAPS 单 stage 运行
- `doc/emseq.md`：EMSeq 入口与配置
- `TEST.md` / `TEST-emseq.md`：维护者回归
