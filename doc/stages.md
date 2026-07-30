# Stages

各 stage 输入输出契约。

TAPS 主流程：

`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> saturation -> summary`

EMSeq 独立入口主流程一致；`mbias` 实现为 `scripts-emseq/mbias.py`（bisulfite 风格，与 TAPS 的 `scripts/mbias.py` 化学规则不同）。详见 `doc/emseq.md`。

## 1. `fastp_split`

**功能**：原始双端 FASTQ 质控与 chunk 切分

**输入**：
- 原始 `R1 FASTQ`
- 原始 `R2 FASTQ`

**输出**：
- `work/<sample>/shard_fastq/*.R1.fq.gz`
- `work/<sample>/shard_fastq/*.R2.fq.gz`
- `work/<sample>/shard_fastq/fastp.html`
- `work/<sample>/shard_fastq/fastp.json`

**约定**：
- 使用 `fastp --split`
- 输出 chunk FASTQ 为 `demux_extract_bc` 直接输入
- 支持 `local` 与 `slurm`
- EMSeq 沿用相同路径与命名

## 2. `demux_extract_bc`

**功能**：从 chunk FASTQ 提取 barcode，分流 matched reads 与 spike-in reads，生成统计

**输入**：
- `shard_fastq/*.R1.fq.gz`
- `shard_fastq/*.R2.fq.gz`

**输出**：
- `demux/<chunk>.R1.demux.fq.gz`
- `demux/<chunk>.R2.demux.fq.gz`
- `demux/<chunk>.R1.spike-in.fq.gz`
- `demux/<chunk>.R2.spike-in.fq.gz`
- `demux/<chunk>.stats.json`

**约定**：
- matched 与 spike-in 必须分流
- read name 改写为 `@barcodeA+barcodeB:original_name`
- 统计文件包含保留率与拒绝原因
- TAPS/EMSeq：`scripts/extract_bc.py`
- 使用 `linker_bc` 定位 barcode：`barcode2(linker_bc左侧)` 与 `barcode1(linker_bc右侧)`（barcode 长度由 whitelist 推导）
- insert 定位：从 `bc1_end` 之后开始搜索 `insert_left`，定位到后对 R1 进行 trim；不再校验中间结构（如 `linker1` / 固定 `others(15 bp)`）

**DBiT-RNA (MVP) 约定**：
- 入口：`scripts-rna/make_cmd.py`
- 实现：`scripts-rna/extract_bc.py`
- `runner=slurm` 时：`#SBATCH --output` / `--error` 默认写入 `work/<sample>/logs/`（本 stage 为 `demux_extract_bc_%x_%j.out` / `.err`）。`align` 为 `align_%x_%j`；`--stage all` 的 driver `run.sbatch` 为 `rna_run_%x_%j`。可在 workflow 的 `slurm.<stage>.output` / `error`、`slurm.run.output` / `error` 或顶层 `slurm_output` / `slurm_error`、CLI `--slurm-output` / `--slurm-error` 覆盖。
- 输入直接使用原始 `r1/r2`（不经过 `fastp_split`）
- RNA 不使用 chunk 概念（单样本单对 clean FASTQ）
- R1 结构：`BC2-linker_bc-BC1-UMI_left-UMI`
- 输出：
  - `demux/<sample>.R1.clean.fq.gz`（序列为 `BC2+BC1+UMI`）
  - `demux/<sample>.R2.clean.fq.gz`（cDNA）
  - `demux/<sample>.stats.json`

## 3. `align`

**功能**：将 demux 后 reads 分别比对至 spike-in 与 host 参考序列

**输入**：
- host：`demux/*.R1.demux.fq.gz`、`demux/*.R2.demux.fq.gz`
- spike-in：`demux/*.R1.spike-in.fq.gz`、`demux/*.R2.spike-in.fq.gz`

**输出**：
- host：`align_shards/<chunk>.cb.bam`
- spike-in：`align_shards/<chunk>.<spike_name>.bam`

**约定**：
- 执行顺序：先 spike-in，后 host
- 支持 0/1/N 个 spike-in
- host 输入必须为 `*.demux.fq.gz`
- spike-in 输入必须为 `*.spike-in.fq.gz`
- EMSeq 使用 `scripts-emseq/aligner.py`（`biscuit align`），输出契约一致

**DBiT-RNA (MVP) 约定**：
- 入口：`scripts-rna/make_cmd.py`
- 实现：`scripts-rna/align.py`
- Slurm 标准输出/错误路径约定见上文 RNA `demux_extract_bc` 小节（本 stage 日志前缀为 `align_`）。
- 输入：`demux/<sample>.R1.clean.fq.gz`、`demux/<sample>.R2.clean.fq.gz`
- 输出：`solo/`（STARsolo matrix 产物目录）
- 附加输出：
  - `solo/star.Solo.out/Gene/raw/barcodes_pos.tsv`
  - `solo/star.Solo.out/Gene/filtered/barcodes_pos.tsv`
  - `barcodes_pos.tsv` 列为 `barcode\tx\ty`，`x/y` 为 whitelist 行号 0-based 坐标
- 当前实现为 STARsolo（matrix-only），不输出 BAM

## 4. `pool`

**功能**：汇总各 chunk 的 host BAM 与 spike-in BAM

**输入**：
- host：`align_shards/*.cb.bam`
- spike-in：`align_shards/*.<spike_name>.bam`

**输出**：
- host：`pooled/pooled.byCB.bam`
- spike-in：`pooled/pooled.<spike_name>.sorted.bam`
- spike-in index：`pooled/pooled.<spike_name>.sorted.bam.bai`

**约定**：
- host 输出固定为 `pooled/pooled.byCB.bam`
- spike-in 输出固定为 `pooled/pooled.<spike_name>.sorted.bam`
- pooled BAM 为 `split`、`mbias`、`call` 上游输入

## 5. `split`

**功能**：按 spot 拆分 pooled host BAM，统计每 spot read 数

**输入**：
- `pooled/pooled.byCB.bam`

**输出**：
- `split_bams/<X_index>/<X_index>_<Y_index>.bam`
- `split_bams/per_spot_read_counts.tsv`

**约定**：
- 按 `CB:Z:<x>+<y>` 解析 spot
- `+` 左侧为 X barcode，右侧为 Y barcode
- `--smoke` 模式最多输出 16 个非空 spot BAM
- `slurm` 模式生成 `commands/05_split_bams.sbatch`，后接排序步骤

## 6. `split` 子步骤 `sort`

**功能**：对 spot BAM 排序与建索引

**输入**：
- `split_bams/**/*.bam`

**输出**：
- `split_bams/**/*.sorted.bam`
- `split_bams/**/*.sorted.bam.bai`

**约定**：
- 为 `split` stage 后处理，不作为顶层 `--stage` 暴露
- 执行内容：`sort + index + remove raw bam`
- `slurm` 模式生成 `commands/05_split_sort.sbatch`
- 单独执行 `split` 时通过 `commands/05_split_submit.sh` 串联依赖

## 7. `mbias`

**功能**：生成 host 和/或 spike-in 的 M-bias 质控结果；为 `call` 准备 host 抽样 BAM

**输入**：
- host：`pooled/pooled.byCB.bam`
- spike-in：`pooled/pooled.<spike_name>.sorted.bam`
- host reference：`call_reference_file`
- spike reference：`spike_in_index`

**输出**：
- `qc/mbias/host.subsampled.sorted.bam`
- `qc/mbias/host.subsampled.sorted.bam.bai`
- `qc/mbias/host.mbias.tsv`
- `qc/mbias/host.mbias.png`
- `qc/mbias/<spike_name>.mbias.tsv`
- `qc/mbias/<spike_name>.mbias.png`

**约定**：
- 默认 `mode=spike`
- 可切换为 `host` 或 `all`
- host 从 `pooled.byCB.bam` 固定比例抽样后排序建索引
- `call` 的 `host_mito` 优先复用 `qc/mbias/host.subsampled.sorted.bam`
- EMSeq：`scripts-emseq/mbias.py` 参考 asTair TOP/BOT 方式，仅统计参考 `CG` 位点；TOP (`99/147`) 以 `C/T` 判甲基化，BOT (`83/163`) 以 `G/A` 判甲基化
- 仅输出 QC，不自动修改 trimming 或 calling 参数
- 输出的 cycle 仅包含 `coverage > 500` 的点，以保证结果稳定性

## 8. `call`

**功能**：生成 host per-spot、host mito 聚合与 spike-in 聚合的甲基化 calling 结果

**输入**：
- host per-spot：`split_bams/**/*.sorted.bam`
- host mito aggregate：优先 `qc/mbias/host.subsampled.sorted.bam`（EMSeq：存在则据此生成 `coverage/host_mito.CG.cov`；否则从 per-spot coverage 汇总线粒体位点）
- spike-in：`pooled/pooled.<spike_name>.sorted.bam`

**输出**：
- `coverage/host/<X_index>/<X_index>_<Y_index>.CG.cov`
- `coverage/host_mito.CG.cov`
- `coverage/<spike_name>.CG.cov`

**约定**：
- `.CG.cov` 仅包含 `coverage>0` 位点
- host 主结果按 spot 输出
- `host_mito` 为单个聚合结果
- `host_mito` 优先复用 `qc/mbias/host.subsampled.sorted.bam`；不存在时从 per-spot `coverage/host/**/*.CG.cov` 汇总线粒体位点（默认 `chrM`）

## 9. `saturation`

**功能**：基于 host per-spot coverage 估计 CpG 饱和度曲线

**输入**：
- `coverage/host/**/*.CG.cov`
- `split_bams/per_spot_read_counts.tsv`

**输出**：
- `qc/saturation/saturation_curve.png`
- `qc/saturation/saturation_summary.tsv`

**约定**：
- 结果写入 `qc/saturation/`
- `summary` 读取 `saturation_summary.tsv` 中的 `saturation_rate`

## 10. `summary`

**功能**：汇总 calling 与 QC 结果为 spot 级与 sample 级报告

**输入**：
- host per-spot calling：`coverage/host/**/*.CG.cov`
- host mito aggregate：`coverage/host_mito.CG.cov`
- spike aggregate calling：`coverage/<spike_name>.CG.cov`
- split reads 统计：`split_bams/per_spot_read_counts.tsv`
- fastp 统计：`shard_fastq/fastp.json`
- demux 统计：`demux/*.stats.json`
- host pooled BAM：`pooled/pooled.byCB.bam`
- spike pooled BAM：`pooled/pooled.<spike_name>.sorted.bam`
- saturation 汇总：`qc/saturation/saturation_summary.tsv`

**输出**：
- `summary/per_spot_summary.tsv`
- `summary/sample_summary.tsv`
- `summary/reads_heatmap.png`
- `summary/cpg_site_count_heatmap.png`
- `summary/mean_methylation_heatmap.png`

**约定**：
- `per_spot_summary.tsv` 每 spot 一行
- 固定列：`X_index`、`Y_index`、`spot`、`mean_methylation`、`cpg_site_count`、`reads`
- `sample_summary.tsv` 每样本一行
- 包含 host、host_mito、各 spike-in 与 `saturation_rate` 汇总
- reads 数值按千分位格式输出
- 缺失输入保持固定列并写 `NA`

## 11. `aggregate`（实验性）

**功能**：将 host per-spot 的 `.CG.cov` 行扁平汇总为两个无表头 TSV（不按坐标跨文件合并；同一 CpG 在不同 spot 文件中仍各占一行）。

**状态**：不纳入 `--stage all`；需显式 `--stage aggregate`。

**输入**：

- `coverage/host/**/*.CG.cov`（与 `call` 产物相同；Bismark-like 六列：`chr`、起始终止坐标、甲基化百分比、`mC`、`C`）

**输出**（均在 `coverage/`，**无表头行**）：

- `aggregated_cg_by_id.tsv`：行序按 `id`（spot，见下）、`chr`、`start`、`end`
- `aggregated_cg_by_pos.tsv`：行序按 `chr`、`start`、`end`、`id`

**列约定**（两文件相同，制表符分隔，顺序固定）：

1. `id`：来源文件名去掉 `.CG.cov` 后缀（例如 `45_21`）
2. `chr`
3. `start`（整数）
4. `end`（整数）
5. `mC`：对应 `.cov` 第 5 列（TAPS `methy_caller` 中为甲基化计数）
6. `C`：对应 `.cov` 第 6 列（非甲基化计数）

**编排**：TAPS：`scripts/make_cmd.py --stage aggregate` 生成本地 `commands/10_aggregate.sh` 或 Slurm `commands/10_aggregate.sbatch`；EMSeq：`scripts-emseq/make_cmd.py --stage aggregate` 生成本地 `commands/11_aggregate.sh` 或 Slurm `commands/11_aggregate.sbatch`（`summary` 已占用 `10_*` 前缀）。等价单步：`scripts/aggregate.py --work-path <work>`。

## 12. `methscan_prepare`（实验性）

**功能**：在独立 pixi 工作区（默认仓库内 `envs/methscan`）中运行 `methscan prepare`，将 host per-spot 的 Bismark-like `*.CG.cov` 转为 methscan 可用的 prepared 输出（仅 host，不含 spike）。

**状态**：不纳入 `--stage all`；需显式 `--stage methscan_prepare`。

**输入**：

- `coverage/host/**/*.CG.cov`（与 `call` 产物相同）

**输出**：

- `methscan/compact/`（由 `methscan prepare` 写入；具体文件以 methscan 版本为准）。规范见 [`docs/developers/contracts.md`](../docs/developers/contracts.md) Methscan 节。

**编排**：TAPS：`scripts/make_cmd.py --stage methscan_prepare` → `commands/11_methscan_prepare.sh` / `11_methscan_prepare.sbatch`；EMSeq：`scripts-emseq/make_cmd.py --stage methscan_prepare` → `commands/12_methscan_prepare.sh` / `12_methscan_prepare.sbatch`。等价单步：`scripts/methscan_run.py prepare --work-path <work>`（可选 `--pixi-manifest <dir>` 覆盖默认 `envs/methscan`；`--chunksize` 默认 `10000000`，与 `methscan prepare` 一致；workflow 键 `methscan_prepare_chunksize` 或 CLI `--methscan-prepare-chunksize` 可覆盖）。
