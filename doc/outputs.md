# Outputs Guide

本页说明 TAPS 主线跑完后先看哪些结果，以及这些结果代表什么。`EMSeq` 独立入口在跑完 `saturation` 与 `summary` 后，同样会产出 `qc/saturation/*` 与 `summary/*`（路径与 TAPS 一致）。

## 推荐检查顺序

1. `summary/sample_summary.tsv`
2. `summary/per_spot_summary.tsv`
3. `summary/*.heatmap.png`
4. `qc/saturation/saturation_summary.tsv`
5. `qc/saturation/saturation_curve.png`

## 主结果清单

### summary

- `work/<sample>/summary/per_spot_summary.tsv`
- `work/<sample>/summary/sample_summary.tsv`
- `work/<sample>/summary/reads_heatmap.png`
- `work/<sample>/summary/cpg_site_count_heatmap.png`
- `work/<sample>/summary/mean_methylation_heatmap.png`

### calling

- `work/<sample>/coverage/host/<X_index>/<X_index>_<Y_index>.CG.cov`
  不包含线粒体 contig（EMSeq 默认去除 `chrM`）
- `work/<sample>/coverage/host_mito.CG.cov`
  为线粒体 contig 的样本级汇总（EMSeq 默认由各 spot host coverage 中的 `chrM` 位点合并得到）
- `work/<sample>/coverage/<spike_name>.CG.cov`

### saturation

- `work/<sample>/qc/saturation/saturation_curve.png`
- `work/<sample>/qc/saturation/saturation_summary.tsv`

### mbias QC

- `work/<sample>/qc/mbias/*.mbias.tsv`
- `work/<sample>/qc/mbias/*.mbias.png`

## 关键字段解释

`per_spot_summary.tsv` 固定包含：

- `X_index`
- `Y_index`
- `spot`
- `mean_methylation`
- `cpg_site_count`
- `reads`

`sample_summary.tsv` 固定列输出样本级结果；缺失输入会写 `NA`，不会改变列结构。除甲基化汇总外，当前还包含：

- `raw_reads`
- `barcoded_reads`
- `barcoded_reads_rate`
- `host_mapped_reads`
- `host_valid_reads`
- `<spike_name>_mapped_reads`
- `valid_reads_rate`
- `saturation_rate`

其中 `saturation_rate` 来自 `qc/saturation/saturation_summary.tsv`；若 `saturation` 未执行或无可用结果，则写 `NA`。

## 常见判断

- `reads` 极低且大面积为 0：优先检查 `demux` 保留率和 `split` 输入
- spike-in 缺失：先检查 `spike_in_index` 配置和 `align/pool` 的 spike 输出
- `host_mito.CG.cov` 缺失：检查 `call` 是否执行，以及 `coverage/host/**/*.CG.cov` 是否已经生成
- `saturation_rate` 为 `NA`：检查 `saturation` 是否执行，以及 `coverage/host/**/*.CG.cov` 和 `split_bams/per_spot_read_counts.tsv` 是否存在

## 相关文档

- `doc/stages.md`：每个 stage 的输入输出契约
- `doc/commands.md`：如何只跑某个 stage
- `TEST.md`：维护者回归检查
