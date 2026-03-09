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
- `workflow`: 测试脚本，用来跑通测试数据

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
python scripts/extract_bc.py \
  "$work_path/shard_fastq/0001.R1.fq.gz" \
  "$work_path/shard_fastq/0001.R2.fq.gz" \
  -b1 configs/barcodes_50a.tsv \
  -b2 configs/barcodes_50a.tsv \
  -o "$work_path/demux/0001"
```

输出：

- `"$work_path/demux/0001.R1.demux.fq.gz"`
- `"$work_path/demux/0001.R1.spike-in.fq.gz"`：未匹配序列
- `"$work_path/demux/0001.R2.demux.fq.gz"`
- `"$work_path/demux/0001.R2.spike-in.fq.gz"`：未匹配序列
- `"$work_path/demux/0001.stats.json"`

匹配到的 reads name 会改为：`@barcodesA+barcodesB:original_reads_name`

3. alignment
4. pool
5. split bam by spots
6. call CpG methylation rates

## 快速开始（命令生成）

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

slurm 的 `fastp_split` 脚本默认包含 `module load fastp`。

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

当 `--runner slurm` 且 `--stage demux_extract_bc` 时，会为每个 chunk 生成一个独立 sbatch 脚本（例如 `02_demux_extract_bc_0001.sbatch`），`--submit` 会逐个提交这些任务。
