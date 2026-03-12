# TEST

本页记录 1.1.2 版本的最小验收和常用回归检查。目标不是覆盖所有内部实现细节，而是快速确认工作流入口、关键 stage 和最终产物是否符合预期。

## 1.1.2 最小验收

建议在提交前至少执行以下 3 步：

1. 检查 CLI 是否可用。
2. 检查 `all` 在 `local` 下能否正确 dry-run。
3. 检查 `all` 在 `slurm` 下能否正确 dry-run。

CLI 检查：

```bash
pixi run python scripts/make_cmd.py --help
```

`local` dry-run：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner local \
  --dry-run
```

`slurm` dry-run：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner slurm \
  --dry-run
```

通过标准：

- 命令正常展开
- stage 顺序为 `fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> summary`
- 无参数缺失、无路径解析错误
- `slurm` dry-run 中外部工具默认应展开为当前 `pixi` 环境下的可执行路径，或显式使用用户传入的 `--*-bin`

## 基础 Smoke

如果只想快速确认入口可用，可以先跑一个最小 smoke：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage demux_extract_bc \
  --runner local \
  --dry-run
```

## 分 Stage 回归

下面这些检查适合在某一 stage 有功能变更时执行。

### `demux_extract_bc`

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

可选：严格基线对照：

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

### `align`

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

可选：验证 spike-in 先于 host：

```bash
pixi run python scripts/align.py \
  --work-path work/test-DNAme-TAPS \
  --chunk 0001 \
  --bwa-index /mnt/e/LiLab_HL/resource/bwa/mm10/genome.fa \
  --spike-in-index lambda=/mnt/e/LiLab_HL/resource/bwa/lambda/genome.fa \
  --spike-in-index puc19=/mnt/e/LiLab_HL/resource/bwa/puc19/genome.fa \
  --dry-run
```

### `pool`

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

### `split`

```bash
pixi run python scripts/split_bams.py --help

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage split \
  --runner local \
  --dry-run

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage split \
  --runner slurm \
  --dry-run
```

可选：`split_bams` smoke：

```bash
pixi run python scripts/split_bams.py \
  --in-bam work/test-DNAme-TAPS/pooled/pooled.byCB.bam \
  --barcodes configs/barcodes_50a.tsv \
  --out-dir work/test-DNAme-TAPS/split_bams.smoke \
  --smoke
```

可选：spot BAM 排序 dry-run：

```bash
pixi run python scripts/bam_sort_parallel.py --help

pixi run python scripts/bam_sort_parallel.py \
  --bam-dir work/test-DNAme-TAPS/split_bams.smoke \
  --jobs 4 \
  --dry-run
```

### `mbias`

```bash
pixi run python scripts/mbias.py --help

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage mbias \
  --runner local \
  --dry-run

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage mbias \
  --runner slurm \
  --dry-run
```

可选：显式分析 host 和 spike-in：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage mbias \
  --mbias-mode all \
  --dry-run
```

检查输出时重点关注：

- `work/<sample>/qc/mbias/*.mbias.tsv`
- `work/<sample>/qc/mbias/*.mbias.png`
- `work/<sample>/qc/mbias/host.subsampled.sorted.bam` 在 `host` 或 `all` 模式下存在
- `lambda` 应接近低甲基化，`puc19` 应保持高甲基化

### `call`

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
```

可选：只检查 host：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage call \
  --call-mode host \
  --dry-run
```

可选：只检查 spike-in：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage call \
  --call-mode spike \
  --dry-run
```

行为检查：

- 若存在 `work/<sample>/qc/mbias/host.subsampled.sorted.bam`，`host_mito` 应直接复用
- 若不存在该 BAM，`call.py` 会先抽样并排序，再生成 `coverage/host_mito.CG.cov`

### `summary`

```bash
pixi run python scripts/summary.py --help

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage summary \
  --dry-run
```

检查输出时重点关注：

- `work/<sample>/summary/per_spot_summary.tsv`
- `work/<sample>/summary/sample_summary.tsv`
- `work/<sample>/summary/reads_heatmap.png`
- `work/<sample>/summary/cpg_site_count_heatmap.png`
- `work/<sample>/summary/mean_methylation_heatmap.png`
- `sample_summary.tsv` 应保持固定列，缺失输入写 `NA`
- 3 张 heatmap 应来自 `per_spot_summary.tsv` 中的 `X_index`、`Y_index` 和对应指标列