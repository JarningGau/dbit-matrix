# Outputs Guide

本页说明跑完后先看哪些结果，以及这些结果代表什么。

## 推荐检查顺序

1. `summary/sample_summary.tsv`  
   先看样本级汇总，确认 host、host_mito、spike-in 是否都产出。

2. `summary/per_spot_summary.tsv`  
   再看 spot 主表，确认每个 spot 的 reads、CpG 位点数和平均甲基化。

3. `summary/*.heatmap.png`  
   最后看空间图，快速识别高低值分布是否符合预期。

## 主结果清单

### summary

- `work/<sample>/summary/per_spot_summary.tsv`
- `work/<sample>/summary/sample_summary.tsv`
- `work/<sample>/summary/reads_heatmap.png`
- `work/<sample>/summary/cpg_site_count_heatmap.png`
- `work/<sample>/summary/mean_methylation_heatmap.png`

### calling

- `work/<sample>/coverage/host/<X_index>/<X_index>_<Y_index>.CG.cov`
- `work/<sample>/coverage/host_mito.CG.cov`
- `work/<sample>/coverage/<spike_name>.CG.cov`

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

`sample_summary.tsv` 固定列输出样本级结果；缺失输入会写 `NA`，不会改变列结构。

## 常见判断

- `reads` 极低且大面积为 0：优先检查 `demux` 保留率和 `split` 输入。
- spike-in 缺失：先检查 `spike_in_index` 配置和 `align/pool` 的 spike 输出。
- `host_mito.CG.cov` 缺失：检查 `call` 是否执行、以及 `mbias` 的 host 抽样 BAM 是否可复用。

## 相关文档

- `doc/stages.md`：想看每个 stage 的输入输出契约
- `doc/commands.md`：想看如何只跑某个 stage
- `TEST.md`：想做维护者回归检查
