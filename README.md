# DBiT-Matrix

该项目旨在为 DBiT-Omics 技术提供统一的上游处理流程，当前覆盖：

- DBiT-RNA
- DBiT-DNAme
  - DBiT-DNAme-TAPS
  - DBiT-DNAme-EMSeq

环境管理使用 `pixi` 和 `nv`，流程需要同时支持 HPC（`sbatch`）与本地 shell。

当前目录说明：

- `scripts`: CLI 脚本
- `data/raw`: 原始数据
- `configs/barcodes_50a.tsv`: spatial barcode whitelist
- `workflow`: 测试脚本，用来跑通测试数据

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

3. alignment
4. pool
5. split bam by spots
6. call CpG methylation rates

## 快速开始（fastp）

1) 使用 workflow 配置生成命令（推荐先 dry-run）

```bash
python scripts/make_cmd.py \
  --workflow-config workflow/fastp_test.json \
  --dry-run
```

不加 `--submit` 时，仅生成命令脚本，不会执行。

2) 生成并立即提交（local）

```bash
python scripts/make_cmd.py \
  --workflow-config workflow/fastp_test.json \
  --submit
```

3) 临时切换到 slurm（CLI 参数优先级高于配置）

```bash
python scripts/make_cmd.py \
  --workflow-config workflow/fastp_test.json \
  --runner slurm \
  --submit
```

slurm 生成脚本默认包含 `module load fastp`。
