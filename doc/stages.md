# Stages

本页描述 dbit-matrix 各个 stage 的输入输出契约。建议把它理解为用户侧的“结果约定”文档：你需要知道每一步读什么、写什么，以及哪些行为是固定的。

固定主流程：

`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> summary`

## 1. `fastp_split`

作用：

- 对原始双端 FASTQ 做基础质控
- 按 chunk 切分，供下游并行处理

输入：

- 原始 `R1 FASTQ`
- 原始 `R2 FASTQ`

输出：

- `work/<sample>/shard_fastq/*.R1.fq.gz`
- `work/<sample>/shard_fastq/*.R2.fq.gz`
- `work/<sample>/fastp.html`
- `work/<sample>/fastp.json`

关键约定：

- 使用 `fastp --split`
- 输出的 chunk FASTQ 是 `demux_extract_bc` 的直接输入
- 支持 `local` 和 `Slurm`

## 2. `demux_extract_bc`

作用：

- 从 chunk FASTQ 中提取 barcode
- 将 matched reads 与 spike-in reads 分流
- 生成保留率与拒绝原因统计

输入：

- `shard_fastq/*.R1.fq.gz`
- `shard_fastq/*.R2.fq.gz`

输出：

- `demux/<chunk>.R1.demux.fq.gz`
- `demux/<chunk>.R2.demux.fq.gz`
- `demux/<chunk>.R1.spike-in.fq.gz`
- `demux/<chunk>.R2.spike-in.fq.gz`
- `demux/<chunk>.stats.json`

关键约定：

- matched 与 spike-in 必须分流
- read name 会被改写为 `@barcodeA+barcodeB:original_name`
- 统计文件必须能反映保留率和拒绝原因

## 3. `align`

作用：

- 将 demux 后的 reads 分别比对到 spike-in 和 host 参考序列

输入：

- host：`demux/*.R1.demux.fq.gz`、`demux/*.R2.demux.fq.gz`
- spike-in：`demux/*.R1.spike-in.fq.gz`、`demux/*.R2.spike-in.fq.gz`

输出：

- host：`align_shards/<chunk>.cb.bam`
- spike-in：`align_shards/<chunk>.<spike_name>.bam`

关键约定：

- 执行顺序固定：先 spike-in，再 host
- 支持 `0/1/N` 个 spike-in
- host 输入必须来自 `*.demux.fq.gz`
- spike-in 输入必须来自 `*.spike-in.fq.gz`

## 4. `pool`

作用：

- 汇总各 chunk 的 host BAM 与 spike-in BAM

输入：

- host：`align_shards/*.cb.bam`
- spike-in：`align_shards/*.<spike_name>.bam`

输出：

- host：`pooled/pooled.byCB.bam`
- spike-in：`pooled/pooled.<spike_name>.sorted.bam`
- spike-in index：`pooled/pooled.<spike_name>.sorted.bam.bai`

关键约定：

- host 最终输出固定为 `pooled/pooled.byCB.bam`
- spike-in 最终输出固定为 `pooled/pooled.<spike_name>.sorted.bam`
- 这些 pooled BAM 是 `split`、`mbias` 和 `call` 的上游输入

## 5. `split`

作用：

- 按 spot 拆分 pooled host BAM
- 统计每个 spot 的 read 数
- 在 `slurm` 编排中，作为 `sort` 的上游步骤单独提交

输入：

- `pooled/pooled.byCB.bam`

输出：

- `split_bams/<X_index>/<X_index>_<Y_index>.bam`
- `split_bams/per_spot_read_counts.tsv`

关键约定：

- 按 `CB:Z:<x>+<y>` 解析 spot
- `+` 左侧是 X barcode，右侧是 Y barcode
- `--smoke` 模式最多输出 16 个非空 spot BAM
- `slurm` 模式生成 `commands/05_split_bams.sbatch`，只负责拆分，不做排序

## 6. `sort`

作用：

- 对 `split` 产生的 spot BAM 做排序和索引
- 在 `slurm` 编排中，必须依赖 `split` 完成后再提交

输入：

- `split_bams/**/*.bam`

输出：

- `split_bams/**/*.sorted.bam`
- `split_bams/**/*.sorted.bam.bai`

关键约定：

- 该步骤是 `split` 的后处理
- 执行内容为 `sort + index + remove raw bam`
- 默认跳过已排序 BAM
- `slurm` 模式生成 `commands/05_split_sort.sbatch`
- 若单独运行 `split` stage，应用 `commands/05_split_submit.sh` 提交，由它负责串联
  `05_split_bams.sbatch -> afterok -> 05_split_sort.sbatch`

## 7. `mbias`

作用：

- 生成 host 和/或 spike-in 的 M-bias 质控结果
- 为 `call` 的聚合 `host_mito` 结果准备可复用的 host 抽样 BAM

输入：

- host：`pooled/pooled.byCB.bam`
- spike-in：`pooled/pooled.<spike_name>.sorted.bam`
- host reference：`call_reference_file`
- spike reference：`spike_in_index`

输出：

- `qc/mbias/host.subsampled.sorted.bam`
- `qc/mbias/host.subsampled.sorted.bam.bai`
- `qc/mbias/host.mbias.tsv`
- `qc/mbias/host.mbias.png`
- `qc/mbias/<spike_name>.mbias.tsv`
- `qc/mbias/<spike_name>.mbias.png`

关键约定：

- 默认 `mode=spike`
- 可显式切换到 `host` 或 `all`
- host 从 `pooled.byCB.bam` 固定比例抽样，再 `sort + index`
- host 抽样 BAM 固定输出为 `qc/mbias/host.subsampled.sorted.bam`
- `call` 的 `host_mito` 会优先复用该 BAM
- spike-in 直接使用全量 `pooled.<spike_name>.sorted.bam`
- 只在参考序列真实 `CpG` 位点上计算甲基化率：`(TG+CA)/(TG+CA+CG)`
- 本步骤只输出 QC，不自动修改 trimming 或 calling 参数

## 8. `call`

作用：

- 生成 host per-spot、host mito 聚合和 spike-in 聚合的甲基化 calling 结果

输入：

- host per-spot：`split_bams/**/*.sorted.bam`
- host mito aggregate：`qc/mbias/host.subsampled.sorted.bam`
- spike-in：`pooled/pooled.<spike_name>.sorted.bam`

输出：

- `coverage/host/<X_index>/<X_index>_<Y_index>.CG.cov`
- `coverage/host_mito.CG.cov`
- `coverage/<spike_name>.CG.cov`

关键约定：

- host 主结果按 spot 输出
- `host_mito` 输出为单个聚合结果
- `host_mito` 优先复用 `qc/mbias/host.subsampled.sorted.bam`
- 若该 BAM 不存在，会从 `pooled/pooled.byCB.bam` 按与 `mbias` 相同规则抽样并排序后再调用

## 9. `summary`

作用：

- 将 calling 结果汇总为 spot 级和 sample 级表格

输入：

- host per-spot calling：`coverage/host/**/*.CG.cov`
- host mito aggregate：`coverage/host_mito.CG.cov`
- spike aggregate calling：`coverage/<spike_name>.CG.cov`
- split reads 统计：`split_bams/per_spot_read_counts.tsv`

输出：

- `summary/per_spot_summary.tsv`
- `summary/sample_summary.tsv`
- `summary/reads_heatmap.png`
- `summary/cpg_site_count_heatmap.png`
- `summary/mean_methylation_heatmap.png`

关键约定：

- `per_spot_summary.tsv` 每个 spot 一行
- 固定包含：`X_index`、`Y_index`、`spot`、`mean_methylation`、`cpg_site_count`、`reads`
- heatmap 固定基于 `per_spot_summary.tsv` 生成
- `reads_heatmap.png` 使用 `(X_index, Y_index, reads)`
- `cpg_site_count_heatmap.png` 使用 `(X_index, Y_index, cpg_site_count)`
- `mean_methylation_heatmap.png` 使用 `(X_index, Y_index, mean_methylation)`
- `sample_summary.tsv` 每个样本一行
- 包含 host、host_mito 和各 spike-in 的汇总甲基化结果
- host spots 平均甲基化按 `cpg_site_count` 加权
- 缺失输入保持固定列并写 `NA`
