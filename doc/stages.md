# Stages

本页描述 `DBiT-Matrix` 各个 stage 的输入输出契约，面向使用者回答“每一步读什么、写什么、哪些行为是固定的”。

TAPS 固定主流程：

`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> saturation -> summary`

EMSeq 当前是独立入口，已覆盖 `fastp_split -> demux_extract_bc -> align -> pool -> split -> call -> saturation -> summary`（不含 `mbias`）。首次使用请看 `doc/emseq.md`。

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
- `work/<sample>/shard_fastq/fastp.html`
- `work/<sample>/shard_fastq/fastp.json`

关键约定：

- 使用 `fastp --split`
- 输出的 chunk FASTQ 是 `demux_extract_bc` 的直接输入
- 支持 `local` 和 `slurm`
- EMSeq 入口沿用相同的输入输出路径和命名

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
- TAPS 使用 `scripts/extract_bc.py`，假定 `R1 = barcodeB-linker2-barcodeA-linker1-Tn5-insert`
- EMSeq 使用 `scripts-emseq/extract_bc.py`，假定 `R1 = linker1-barcodeB-linker2-barcodeA-others(15 bp)-Tn5-insert`
- EMSeq demux 的外部参数设计与 TAPS 对齐；`others(15 bp)` 作为内部固定约束，不单独暴露为 CLI 参数

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
- EMSeq 入口使用 `scripts-emseq/aligner.py`，以 `biscuit align` 产出比对结果，并保持相同的输出契约：host 为 `align_shards/<chunk>.cb.bam`，spike-in 为 `align_shards/<chunk>.<spike_name>.bam`

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

输入：

- `pooled/pooled.byCB.bam`

输出：

- `split_bams/<X_index>/<X_index>_<Y_index>.bam`
- `split_bams/per_spot_read_counts.tsv`

关键约定：

- 按 `CB:Z:<x>+<y>` 解析 spot
- `+` 左侧是 X barcode，右侧是 Y barcode
- `--smoke` 模式最多输出 16 个非空 spot BAM
- `slurm` 模式先生成 `commands/05_split_bams.sbatch`，再串联排序步骤

## 6. `split` 内部子步骤 `sort`

作用：

- 对 `split` 产生的 spot BAM 做排序和索引

输入：

- `split_bams/**/*.bam`

输出：

- `split_bams/**/*.sorted.bam`
- `split_bams/**/*.sorted.bam.bai`

关键约定：

- 该步骤是 `split` stage 的后处理，不作为顶层 `--stage` 暴露
- 执行内容为 `sort + index + remove raw bam`
- `slurm` 模式生成 `commands/05_split_sort.sbatch`
- 单独运行 `split` 时可通过 `commands/05_split_submit.sh` 串联依赖

## 7. `mbias`

作用：

- 生成 host 和或 spike-in 的 M-bias 质控结果
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
- host 从 `pooled.byCB.bam` 固定比例抽样，再排序并建立索引
- `call` 的 `host_mito` 会优先复用 `qc/mbias/host.subsampled.sorted.bam`
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

- `.CG.cov` 仅包含 `coverage>0` 的位点行
- host 主结果按 spot 输出
- `host_mito` 输出为单个聚合结果
- `host_mito` 优先复用 `qc/mbias/host.subsampled.sorted.bam`

## 9. `saturation`

作用：

- 基于 host per-spot coverage 估计样本的 CpG 饱和度曲线

输入：

- `coverage/host/**/*.CG.cov`
- `split_bams/per_spot_read_counts.tsv`

输出：

- `qc/saturation/saturation_curve.png`
- `qc/saturation/saturation_summary.tsv`

关键约定：

- 结果写到 `qc/saturation/`
- `summary` 会读取 `saturation_summary.tsv` 中的 `saturation_rate`

## 10. `summary`

作用：

- 将 calling 和 QC 结果汇总为 spot 级和 sample 级结果

输入：

- host per-spot calling：`coverage/host/**/*.CG.cov`
- host mito aggregate：`coverage/host_mito.CG.cov`
- spike aggregate calling：`coverage/<spike_name>.CG.cov`
- split reads 统计：`split_bams/per_spot_read_counts.tsv`
- fastp 统计：`shard_fastq/fastp.json`
- demux 统计：`demux/*.stats.json`
- host pooled BAM：`pooled/pooled.byCB.bam`
- spike pooled BAM：`pooled/pooled.<spike_name>.sorted.bam`
- saturation 汇总：`qc/saturation/saturation_summary.tsv`

输出：

- `summary/per_spot_summary.tsv`
- `summary/sample_summary.tsv`
- `summary/reads_heatmap.png`
- `summary/cpg_site_count_heatmap.png`
- `summary/mean_methylation_heatmap.png`

关键约定：

- `per_spot_summary.tsv` 每个 spot 一行
- 固定包含：`X_index`、`Y_index`、`spot`、`mean_methylation`、`cpg_site_count`、`reads`
- `sample_summary.tsv` 每个样本一行
- 包含 host、host_mito、各 spike-in 和 `saturation_rate` 的汇总结果
- reads 数值字段按千分位格式输出
- 缺失输入保持固定列并写 `NA`
