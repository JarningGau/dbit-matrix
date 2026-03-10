# DBiT-Matrix

该项目旨在为 DBiT-Omics 技术提供统一的上游处理流程，当前覆盖：

- DBiT-RNA
- DBiT-DNAme
  - DBiT-DNAme-TAPS
  - DBiT-DNAme-EMSeq

环境管理使用 `pixi + uv`，流程需要同时支持 HPC（`sbatch`）与本地 shell。

当前目录说明：

- `scripts`: CLI 脚本
- `data/raw`: 原始数据
- `configs/barcodes_50a.tsv`: spatial barcode whitelist
- `workflow/dbit_taps_test.json`: 参数配置文件，参数文件命名方式为`<方法>_<样本名>`

## 环境固定（pixi + uv）

首次初始化并生成锁文件：

```bash
pixi lock
pixi install
```

之后运行脚本统一使用 `pixi run`：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --dry-run
```

如果未来需要增加 PyPI 依赖（由 `uv` 解析并写入锁文件），使用：

```bash
pixi add --pypi <package>
pixi lock
```

## 环境测试（验收）

建议在每次拉取更新后执行以下检查，确认环境可用且可复现。

1) 安装环境

```bash
pixi install
```

2) 检查关键工具来源和版本

```bash
pixi run which python
pixi run python --version
pixi run which fastp
pixi run fastp --version
pixi run which bwa
pixi run bwa 2>&1 | head -n 1
pixi run which sinto
pixi run sinto --help >/dev/null
pixi run which samtools
pixi run samtools --version | head -n 1
```

3) 运行最小 smoke test（不执行实际提交）

```bash
pixi run make-cmd-dry-run
```

4) 验证 lockfile 可复现

在另一台机器（或清理本地 `.pixi` 后）仅保留 `pixi.toml + pixi.lock`，重新执行：

```bash
pixi install
pixi run make-cmd-dry-run
```

若命令可正常执行且输出一致，可认为环境固定成功。

## DBiT-DNAme-TAPS

测试数据：

- R1: `data/raw/test-DNAme-TAPS/250812_intestine_dbit_taps_R1.fq.gz`
- R2: `data/raw/test-DNAme-TAPS/250812_intestine_dbit_taps_R2.fq.gz`

文库结构：

- R1: `barcodeB-linker2-barcodeA-linker1-Tn5-insert`
- R2: `insert`

参数：

- spatial barcodes (barcodeA and barcodeB): `configs/barcodes_50a.tsv`
- linker1: `GTGGCCGATGTTTCG`
- linker2: `ATCCACGTGCTTGAGAGGCCAGAGCATTCG`
- tn5: `CATCGGCGTACGACTAGATGTGTATAAGAGACAG`
- linker_edit_distance: linker1/linker2/tn5 的模糊匹配编辑距离上限（默认 `1`）
- barcode_hamming_distance: barcode whitelist 纠错的汉明距离上限（默认 `1`）
- gzip_level: demux 输出 FASTQ 的 gzip 压缩等级（`0-9`，默认 `6`）
- bwa_index: `bwa mem` 使用的索引前缀（align 阶段必填）
- bwa_threads: align 阶段 `bwa mem -t` 线程数（默认 `2`）
- bwa_bin / sinto_bin: align 阶段可执行程序名或路径（默认 `bwa` / `sinto`）
- samtools_bin: spike-in 对齐输出 BAM 使用的 `samtools` 程序（默认 `samtools`）
- spike_in_index: spike-in 参考配置，可写成对象（如 `{ "lambda": "...", "puc19": "..." }`）或 `["lambda=...", "puc19=..."]`

## 核心处理流程

1. split FASTQ to chunks

```bash
fastp \
  -i "$R1" \
  -I "$R2" \
  -o "$work_path/shard_fastq/R1.fq.gz" \
  -O "$work_path/shard_fastq/R2.fq.gz" \
  -w "$fastp_threads" \
  --split "$number_of_split_parts" \
  --disable_adapter_trimming \
  -h "$work_path/shard_fastq/fastp.html" \
  -j "$work_path/shard_fastq/fastp.json"
```

2. demultiplexing (demux): extract barcode from R1  
   每个 chunk 一个命令：

```bash
pixi run python scripts/extract_bc.py \
  "$work_path/shard_fastq/0001.R1.fq.gz" \
  "$work_path/shard_fastq/0001.R2.fq.gz" \
  -b1 configs/barcodes_50a.tsv \
  -b2 configs/barcodes_50a.tsv \
  -o "$work_path/demux/0001" \
  --linker-edit-distance 1 \
  --barcode-hamming-distance 1
```

输出：

- `"$work_path/demux/0001.R1.demux.fq.gz"`
- `"$work_path/demux/0001.R1.spike-in.fq.gz"`：未匹配序列
- `"$work_path/demux/0001.R2.demux.fq.gz"`
- `"$work_path/demux/0001.R2.spike-in.fq.gz"`：未匹配序列
- `"$work_path/demux/0001.stats.json"`

匹配到的 reads name 会改为：`@barcodesA+barcodesB:original_reads_name`
运行 `extract_bc.py` 时会周期性输出处理速度（`xx reads/s`），默认每秒刷新一次；
可用 `--progress-interval-seconds` 调整刷新间隔（设为 `0` 可关闭）。
可用 `--gzip-level` 调整输出压缩等级（`0-9`，等级越低通常越快）。
定位和匹配策略：

- 模糊匹配底层使用 `fuzzysearch.find_near_matches`（仅在 exact 失败时启用）。
- `linker2` 在前端窗口内定位（`0 ~ bc2_len + len(linker2) + 12`），并要求窗口内唯一可接受命中（多命中或并列最优视为歧义，直接拒绝）。
- 用 linker2 左侧/右侧固定长度提取 barcodeB/barcodeA，再校验 linker1。
- `linker1` 采用分层容错：exact -> mismatch-only -> full edit distance。
- Tn5 检索用于界定 insert 起点：
  - `len(tn5) > 15` 时仅用 Tn5 右侧 15bp 检索；
  - 否则使用 Tn5 全长检索；
  - 二者都只在 `after_pos` 后的小窗口内搜索并要求唯一可接受命中，命中后直接计算 `tn5_end`。
- whitelist 匹配顺序为：精确匹配 -> 汉明距离退化匹配（需唯一最优）。

3. alignment

software: `bwa`, `sinto`, `samtools`

spike-in 可以有多个参考基因组（例如 lambda、pUC19），
建议每个参考单独比对并按 `<chunk>.<spike_name>.bam` 命名输出。
在 `scripts/align.py` 中执行顺序固定为：先对每个 spike-in 参考比对，再对 host genome 比对并生成 `.cb.bam`。

```bash
bwa mem -t 2 \
$bwa_index_host_genome \
$work_path/demux/0001.R1.demux.fq.gz \
$work_path/demux/0001.R2.demux.fq.gz |\
sinto nametotag -b - -O b -o $work_path/align_shards/0001.cb.bam

bwa mem -t 2 \
$bwa_index_lambda \
$work_path/demux/0001.R1.spike-in.fq.gz \
$work_path/demux/0001.R2.spike-in.fq.gz |\
samtools view -b -o $work_path/align_shards/0001.lambda.bam -

bwa mem -t 2 \
$bwa_index_puc19 \
$work_path/demux/0001.R1.spike-in.fq.gz \
$work_path/demux/0001.R2.spike-in.fq.gz |\
samtools view -b -o $work_path/align_shards/0001.puc19.bam -
```

4. pool

输入来源：
- host: `"$work_path/align_shards/"*.cb.bam`
- spike-in: `"$work_path/align_shards/"*.<spike_name>.bam`

执行约定：
- `runner=local`：一个脚本内先处理 spike-in，再处理 host。
- `runner=slurm`：拆成两个 sbatch，分别处理 spike-in 与 host。

```bash
## host genome
samtools cat -o $work_path/pooled/pooled.cb.bam $work_path/align_shards/*.cb.bam
samtools sort -m 16G -@ 4 -t CB -o $work_path/pooled/pooled.byCB.bam $work_path/pooled/pooled.cb.bam
rm $work_path/pooled/pooled.cb.bam

## spike-in
## pool -> sort by coords -> index -> remove tmp file
### lambda
samtools cat -o $work_path/pooled/pooled.lambda.bam $work_path/align_shards/*.lambda.bam && samtools sort -o $work_path/pooled/pooled.lambda.sorted.bam $work_path/pooled/pooled.lambda.bam && samtools index $work_path/pooled/pooled.lambda.sorted.bam && rm $work_path/pooled/pooled.lambda.bam
### puc19
samtools cat -o $work_path/pooled/pooled.puc19.bam $work_path/align_shards/*.puc19.bam && samtools sort -o $work_path/pooled/pooled.puc19.sorted.bam $work_path/pooled/pooled.puc19.bam && samtools index $work_path/pooled/pooled.puc19.sorted.bam && rm $work_path/pooled/pooled.puc19.bam
```

5. split bam by spots
6. call CpG methylation rates

## 生成 fastp 命令

1) 使用 workflow 配置生成命令（推荐先 dry-run）

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --dry-run
```

不加 `--submit` 时，仅生成命令脚本，不会执行。

2) 生成并立即提交（local）

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --submit
```

3) 临时切换到 slurm（CLI 参数优先级高于配置）

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --runner slurm \
  --submit
```

slurm 脚本使用 `pixi` 环境中的工具（不依赖 `module load`）。

## 生成 demux 命令

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage demux_extract_bc \
  --dry-run
```

生成并立刻提交
```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage demux_extract_bc \
  --submit
```

执行时会扫描 `"$work_path/shard_fastq/"*.R1.fq.gz`，并按 chunk 调用 `scripts/extract_bc.py`。
当 `--runner local` 时，demux 执行会显示按 chunk 的进度条。

当 `--runner slurm` 且 `--stage demux_extract_bc` 时，会为每个 chunk 生成一个独立 sbatch 脚本（例如 `02_demux_extract_bc_0001.sbatch`），`--submit` 会逐个提交这些任务。

## 生成 align 命令

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage align \
  --dry-run
```

生成并立刻提交
```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage align \
  --submit
```

执行时会扫描 `"$work_path/demux/"*.R1.demux.fq.gz`，并要求存在对应的 `*.R2.demux.fq.gz`。
当 `--runner local` 时，生成 `03_align.sh`，内部按 chunk 顺序运行 `scripts/align.py`；每个 chunk 会在同一条执行链中先跑全部 spike-in，再跑 host。
当 `--runner slurm` 且 `--stage align` 时，会为每个 chunk 生成一个独立 sbatch 脚本（例如 `03_align_0001.sbatch`）；单个 sbatch 内绑定该 chunk 的 spike-in + host 对齐顺序执行，`--submit` 会逐个提交这些任务。

## 生成 pool 命令

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage pool \
  --dry-run
```

生成并立刻提交
```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage pool \
  --submit
```

执行时会扫描 `"$work_path/align_shards/"*.cb.bam` 和 `"$work_path/align_shards/"*.<spike_name>.bam`。  
当 `--runner local` 时，生成 `04_pool.sh`，内部先处理 spike-in，再处理 host。  
当 `--runner slurm` 且 `--stage pool` 时，会生成两个 sbatch：`04_pool_spike.sbatch` 与 `04_pool_host.sbatch`；`--submit` 会依次提交这两个任务。
