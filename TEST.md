# TEST (Maintainer)

本页是维护者/贡献者回归文档，不是首次用户上手教程。  
目标是快速确认工作流入口、关键 stage 和最终产物在版本演进后仍符合预期。

如果你是首次使用者，请先阅读：`README.md -> doc/setup.md -> doc/config.md -> doc/commands.md -> doc/outputs.md`。

## 最小验收

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
- stage 顺序为 `fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> saturation -> summary`
- 无参数缺失、无路径解析错误
- `slurm` dry-run 中外部工具默认应展开为当前 `pixi` 环境下的可执行路径，或显式使用用户传入的 `--*-bin`

`slurm` 真实生成回归：

```bash
rm -rf work

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner slurm
```

通过标准：

- 空 `work/` 下也能成功生成，不要求预先存在 `shard_fastq`、`demux`、`split_bams` 等下游输入
- 会写出 `work/<sample>/commands/run.sbatch`
- 会直接写出各 stage/chunk 的 sbatch（例如 `02_demux_extract_bc_0001.sbatch`、`03_align_0001.sbatch`），不依赖 `run_02_*.sbatch` launcher 链

`slurm` 一次提交回归（无 nested submit）：

```bash
rm -rf work

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner slurm \
  --submit
```

通过标准：

- 提交流程由客户端一次完成（日志包含 `submit_mode=client_side_sbatch_dag`）
- 不会创建 `dbit_all_launcher_*` 类型的 launcher 作业

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

检查点：

- `slurm` 模式下应显示两个 stage 脚本：`commands/05_split_bams.sbatch`、`commands/05_split_sort.sbatch`
- 真实生成时应额外写出 `commands/05_split_submit.sh`，并由它负责以 `afterok` 串联两个 `sbatch`

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
- trimming 应按原始 read 两端 cycle 生效，而不是按 BAM `query_pos` 固定方向裁切；reverse-strand reads 的 left/right 语义应与 `mbias` 一致
- `scripts/methy_caller.py` 应按 batch 流式写出 `.CG.cov`，不再先累计完整结果再统一落盘
- `.CG.cov` 文件不应包含 `coverage=0` 的行

### `saturation`

```bash
pixi run python scripts/saturation.py --help

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage saturation \
  --runner local \
  --dry-run

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage saturation \
  --runner slurm \
  --dry-run
```

检查输出时重点关注：

- `work/<sample>/qc/saturation/saturation_curve.png`
- `work/<sample>/qc/saturation/saturation_summary.tsv`
- `saturation_summary.tsv` 应包含：`sample_id`、`observed_median_unique_cpgs`、
  `theoretical_max_median_unique_cpgs`、`predicted_median_unique_cpgs_at_2x`、
  `saturation_rate`、`hq_spot_count`

### `summary`

```bash
pixi run python scripts/summary.py --help

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage summary \
  --dry-run

pixi run python scripts/summary.py \
  --work-path work/test-DNAme-TAPS \
  --spike-in-name lambda \
  --spike-in-name puc19 \
  --dry-run
```

检查输出时重点关注：

- `work/<sample>/summary/per_spot_summary.tsv`
- `work/<sample>/summary/sample_summary.tsv`
- `work/<sample>/summary/reads_heatmap.png`
- `work/<sample>/summary/cpg_site_count_heatmap.png`
- `work/<sample>/summary/mean_methylation_heatmap.png`
- `sample_summary.tsv` 应保持固定列，缺失输入写 `NA`
- `sample_summary.tsv` 应包含：`raw_reads`、`barcoded_reads`、`host_mapped_reads`、
  `lambda_mapped_reads`、`puc19_mapped_reads`、`host_valid_reads`、
  `barcoded_reads_rate`、`valid_reads_rate`（若输入缺失则为 `NA`）
- `barcoded_reads` 应按 reads 口径（demux 的 read pairs 计数需要 `×2` 归一化）
- reads 数值字段应为千分位格式，百分比字段应为 `XX.XX%` 格式
- `valid_reads_rate` 应满足 `host_valid_reads/raw_reads`
- 3 张 heatmap 应来自 `per_spot_summary.tsv` 中的 `X_index`、`Y_index` 和对应指标列