# Stages

## 1. `fastp_split`

输入：

- 原始 `R1/R2 FASTQ`

输出：

- `work/<sample>/shard_fastq/*.R1.fq.gz`
- `work/<sample>/shard_fastq/*.R2.fq.gz`
- `fastp.html`
- `fastp.json`

说明：

- 用 `fastp --split` 将原始 FASTQ 切成 chunk
- 当前支持本地与 Slurm

## 2. `demux_extract_bc`

输入：

- `shard_fastq/*.R1.fq.gz`
- `shard_fastq/*.R2.fq.gz`

输出：

- `demux/<chunk>.R1.demux.fq.gz`
- `demux/<chunk>.R2.demux.fq.gz`
- `demux/<chunk>.R1.spike-in.fq.gz`
- `demux/<chunk>.R2.spike-in.fq.gz`
- `demux/<chunk>.stats.json`

说明：

- matched 与 spike-in 必须分流
- reads name 会改写为 `@barcodeA+barcodeB:original_name`
- stats 需反映保留率与拒绝原因

## 3. `align`

输入：

- host: `demux/*.R1/2.demux.fq.gz`
- spike-in: `demux/*.R1/2.spike-in.fq.gz`

输出：

- host: `align_shards/<chunk>.cb.bam`
- spike-in: `align_shards/<chunk>.<spike_name>.bam`

说明：

- 执行顺序固定：先 spike-in，再 host
- 支持 0/1/N 个 spike-in
- Slurm 下每个 chunk 一个 sbatch

## 4. `pool`

输入：

- host: `align_shards/*.cb.bam`
- spike-in: `align_shards/*.<spike_name>.bam`

输出：

- host: `pooled/pooled.byCB.bam`
- spike-in: `pooled/pooled.<spike_name>.sorted.bam` 及索引

说明：

- 本地模式：先 spike-in，再 host
- Slurm 模式：拆成 `spike` 与 `host` 两个 job

## 5. `split`

输入：

- `pooled/pooled.byCB.bam`

输出：

- `split_bams/<X_index>/<X_index>_<Y_index>.bam`
- `split_bams/per_spot_read_counts.tsv`

说明：

- 当前 `CB` 格式：`CB:Z:<x>+<y>`
- `+` 左侧为 X barcode，右侧为 Y barcode
- `--smoke` 仅随机输出最多 16 个非空 spot BAM

## 6. `sort`（split 后处理）

输入：

- `split_bams/**/*.bam`

输出：

- `split_bams/**/*.sorted.bam`
- `split_bams/**/*.sorted.bam.bai`

说明：

- 执行 `sort + index + remove raw bam`
- 默认跳过已排序 BAM
- Slurm 下与 `split_bams` 分成不同资源配置

## 7. `mbias`

状态：

- 已实现（MVP）

输入：

- host: `pooled/pooled.byCB.bam`
- spike-in: `pooled/pooled.<spike_name>.sorted.bam`

输出：

- `qc/mbias/host.subsampled.sorted.bam`
- `qc/mbias/host.subsampled.sorted.bam.bai`
- `qc/mbias/host.mbias.tsv`
- `qc/mbias/<spike_name>.mbias.tsv`

说明：

- host 从 `pooled.byCB.bam` 固定比例抽样（固定 seed），随后 `sort + index` 再做 mbias
- spike-in 直接使用全量 `pooled.<spike_name>.sorted.bam` 做 mbias
- 仅输出 QC 结果，不自动修改 trimming/calling 参数

## 8. `call`

状态：

- 已实现（MVP）

输入：

- host by spots: `split_bams/**/*.sorted.bam`
- spike-in: `pooled/pooled.<spike_name>.sorted.bam`

输出：

- `coverage/host/<X_index>/<X_index>_<Y_index>.CG.cov`
- `coverage/host_mito/<X_index>/<X_index>_<Y_index>.CG.cov`
- `coverage/lambda.CG.cov`
- `coverage/puc19.CG.cov`

说明：

- host 按 spot 输出；每个 spot 额外输出一份线粒体位点结果
- 本地模式：单脚本内并行 spots，spike-in 顺序执行
- Slurm 模式：host 只生成一个 sbatch，在作业内并行处理 spots；spike-in 每个 reference 一个 sbatch
