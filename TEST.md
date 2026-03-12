# TEST

## 当前测试覆盖

- 已覆盖阶段：`demux_extract_bc`、`align`、`pool`、`split`、`mbias`、`call`
- 当前未覆盖完整实现：`split -> mbias -> call` 端到端验收
- 建议提交前至少执行：CLI `--help` 检查 + `scripts/make_cmd.py --dry-run` 回归

## 基础 smoke

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage demux_extract_bc \
  --runner local \
  --dry-run
```

## Demux 回归

最小功能回归：

```bash
pixi run python scripts/extract_bc.py \
  work/test-DNAme-TAPS/shard_fastq/0001.R1.fq.gz \
  work/test-DNAme-TAPS/shard_fastq/0001.R2.fq.gz \
  -b1 configs/barcodes_50a.tsv \
  -b2 configs/barcodes_50a.tsv \
  -o work/test-DNAme-TAPS/demux/0001.test \
  --linker-edit-distance 1 \
  --barcode-hamming-distance 1 \
  --gzip-level 1
```

参考输出摘要：

```text
[extract_bc] kept=346137/487575 spike_in=141438 avg_speed=41021.6 reads/s
```

可选：严格基线（仅精确匹配）用于速度/结果对照：

```bash
pixi run python scripts/extract_bc.py \
  work/test-DNAme-TAPS/shard_fastq/0001.R1.fq.gz \
  work/test-DNAme-TAPS/shard_fastq/0001.R2.fq.gz \
  -b1 configs/barcodes_50a.tsv \
  -b2 configs/barcodes_50a.tsv \
  -o work/test-DNAme-TAPS/demux/0001.strict \
  --linker-edit-distance 0 \
  --barcode-hamming-distance 0 \
  --gzip-level 1
```

参考输出摘要：

```text
[extract_bc] kept=275554/487575 spike_in=212021 avg_speed=64884.6 reads/s
```

## Align 回归

CLI 帮助与命令生成：

```bash
pixi run python scripts/align.py --help

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage align \
  --runner local \
  --dry-run

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage align \
  --runner slurm \
  --dry-run
```

可选：验证多 spike-in 顺序（先 spike-in 再 host）：

```bash
pixi run python scripts/align.py \
  --work-path work/test-DNAme-TAPS \
  --chunk 0001 \
  --bwa-index /mnt/e/LiLab_HL/resource/bwa/mm10/genome.fa \
  --spike-in-index lambda=/mnt/e/LiLab_HL/resource/bwa/lambda/genome.fa \
  --spike-in-index puc19=/mnt/e/LiLab_HL/resource/bwa/puc19/genome.fa \
  --dry-run
```

## Pool 回归

```bash
pixi run python scripts/pool.py --help

pixi run python scripts/pool.py \
  --work-path work/test-DNAme-TAPS \
  --samtools-bin samtools \
  --spike-in-name lambda \
  --spike-in-name puc19 \
  --dry-run

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage pool \
  --runner local \
  --dry-run

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage pool \
  --runner slurm \
  --dry-run
```

## Split 回归

```bash
pixi run python scripts/split_bams.py --help

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage split \
  --dry-run

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage split \
  --runner slurm \
  --dry-run
```

可选：`split_bams` smoke 模式，仅随机输出 16 个非空 spot BAM：

```bash
pixi run python scripts/split_bams.py \
  --in-bam work/test-DNAme-TAPS/pooled/pooled.byCB.bam \
  --barcodes configs/barcodes_50a.tsv \
  --out-dir work/test-DNAme-TAPS/split_bams.smoke \
  --smoke
```

可选：spot BAM 并行排序 dry-run：

```bash
pixi run python scripts/bam_sort_parallel.py --help

pixi run python scripts/bam_sort_parallel.py \
  --bam-dir work/test-DNAme-TAPS/split_bams.smoke \
  --jobs 4 \
  --dry-run
```

## Call 回归

```bash
pixi run python scripts/methy_caller.py --help

pixi run python scripts/call.py --help

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage call \
  --call-mode all \
  --call-r1-left-trimming 0 \
  --call-r1-right-trimming 0 \
  --call-r2-left-trimming 0 \
  --call-r2-right-trimming 0 \
  --dry-run

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage call \
  --runner slurm \
  --call-mode all \
  --call-r1-left-trimming 0 \
  --call-r1-right-trimming 0 \
  --call-r2-left-trimming 0 \
  --call-r2-right-trimming 0 \
  --dry-run

# 可选：只生成 host calling 命令
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage call \
  --call-mode host \
  --dry-run

# 可选：只生成 spike calling 命令
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage call \
  --call-mode spike \
  --dry-run
```

## Mbias 回归

```bash
pixi run python scripts/mbias.py --help

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage mbias \
  --dry-run

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage mbias \
  --runner slurm \
  --dry-run

# 可选：显式分析 host + spike
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage mbias \
  --mbias-mode all \
  --dry-run
```

检查输出（非 dry-run）：

- `work/<sample>/qc/mbias/*.mbias.tsv`
- `work/<sample>/qc/mbias/*.mbias.png`
- `mbias` 结果应与 `call` 方向一致：`lambda` 接近低甲基化，`puc19` 保持高甲基化

可选：验证 `call.py` 主机位点并行 dry-run：

```bash
pixi run python scripts/call.py \
  --work-path work/test-DNAme-TAPS \
  --mode host \
  --reference-file /mnt/e/LiLab_HL/resource/bwa/mm10/genome.fa \
  --chromosomes chr1,chr2,chr3 \
  --mito-chromosomes chrM \
  --jobs 4 \
  --r1-left-trimming 5 \
  --r1-right-trimming 5 \
  --r2-left-trimming 5 \
  --r2-right-trimming 5 \
  --dry-run
```

## 下一里程碑

- 在 `split -> mbias -> call` 之间补齐端到端验收命令