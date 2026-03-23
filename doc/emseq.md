# EMSeq Guide

本页是 `EMSeq` 独立入口的用户说明，覆盖当前 MVP 范围。维护者回归检查请看 `TEST-emseq.md`。

## 当前支持范围

- 当前状态：试验性 / MVP
- 入口脚本：`scripts-emseq/make_cmd.py`
- 当前支持 stage：`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> saturation -> summary`，以及 workflow stage **`all`**（生成 `commands/run.sh` 或 `commands/run.sbatch`，按该顺序串联各 stage；与下面 `mbias_mode` / `call_mode` 的 `all` 含义不同）
- `mbias` 单步脚本：`scripts-emseq/mbias.py`（bisulfite 风格 M-bias；参考 asTair 的 TOP/BOT 方式，仅统计 reference `CG` 位点；与 TAPS 的 `scripts/mbias.py` 判定规则不同）

## 与 TAPS 的主要差异

- `EMSeq` 使用独立入口，不通过 `scripts/make_cmd.py`
- `demux_extract_bc` 使用 `scripts-emseq/extract_bc.py`
- `align` 使用 `scripts-emseq/aligner.py`
- EMSeq 假定的 `R1` 结构为：
  `linker1-barcodeB-linker2-barcodeA-others(15 bp)-Tn5-insert`
- `align` 使用 `biscuit align` 进行比对，并将 host 结果整理为带 `CB` tag 的 BAM
- demux 和 align 的外部参数风格与 TAPS 对齐，当前仍通过 `workflow/*.json` 配置

## 最小配置

建议从 `workflow/dbit_emseq_test.json` 复制，并至少确认以下字段：

- `sample_id`
- `r1`
- `r2`
- `work_root`
- `fastp_threads`
- `number_of_split_parts`
- `fastp_bin`
- `barcode1_whitelist`
- `barcode2_whitelist`
- `linker1`
- `linker2`
- `tn5`
- `biscuit_reference`
- `biscuit_threads`
- `biscuit_batch_size`
- `biscuit_bin`
- `sinto_bin`
- `samtools_bin`

可调的 demux 参数：

- `linker_edit_distance`
- `barcode_hamming_distance`
- `gzip_level`

可调的 align 参数：

- `spike_in_index`

可调的 mbias 参数：

- `mbias_script`：默认 `scripts-emseq/mbias.py`
- `mbias_mode`：`all` / `host` / `spike`（默认 `spike`）；`all` 或 `host` 需要 `call_reference_file`；`all` 或 `spike` 需要 `spike_in_index`
- `mbias_host_subsample_fraction`、`mbias_max_cycle`、`mbias_min_mapping_quality`
- Slurm：`slurm.mbias.host` / `slurm.mbias.spike`（partition / mem / cpus_per_task）

可调的 call 参数：

- `call_mito_chromosomes`：host 中哪些 contig 视为线粒体。默认 `chrM`；会从 per-spot `coverage/host/**/*.CG.cov` 中剔除这些 contig；`coverage/host_mito.CG.cov` 优先由 `qc/mbias/host.subsampled.sorted.bam` 经 pileup 得到（若该 BAM 不存在，则从 per-spot 文件中汇总的线粒体位点生成）。

可调的 saturation 参数（复用 `scripts/saturation.py`，在 `call` 之后运行）：

- `saturation_script`：默认 `scripts/saturation.py`
- `saturation_reads_threshold`：HQ spot 的 reads 阈值，默认 `1e6`
- Slurm：`slurm.saturation`（partition / mem / cpus_per_task）

可调的 summary 参数（复用 `scripts/summary.py`，在 `saturation` 之后运行）：

- `summary_script`：默认 `scripts/summary.py`
- `spike_in_index`：用于向 `summary` 传递 `--spike-in-name`（与 `spike_in_index` 的 key 一致）；未配置时由 `summary` 从 `coverage/` 下发现 spike 覆盖文件
- Slurm：`slurm.summary`（partition / mem / cpus_per_task）

## 首次运行命令

一键检查整条主线（推荐先 `--dry-run`，不落盘）：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage all \
  --runner local \
  --dry-run
```

Slurm 下同样可先 dry-run，再生成 `commands/run.sbatch` 驱动（`--runner slurm` 且非 dry-run 时写入；`all` 驱动脚本的 SBATCH 资源默认取自 `slurm.summary`）：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage all \
  --runner slurm \
  --dry-run
```

先检查 `fastp_split`：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage fastp_split \
  --runner local \
  --dry-run
```

再检查 `demux_extract_bc`：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage demux_extract_bc \
  --runner local \
  --dry-run
```

如需提交到 Slurm，把 `--runner local` 改为 `--runner slurm`。

最后检查 `align`：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage align \
  --runner local \
  --dry-run
```

再检查 `pool`：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage pool \
  --runner local \
  --dry-run
```

再检查 `split`：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage split \
  --runner local \
  --dry-run
```

再检查 `mbias`（需已有 `pooled/` BAM；`mbias_mode` 为 `host` 时需配置 `call_reference_file`）：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage mbias \
  --runner local \
  --dry-run
```

最后检查 `call`：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage call \
  --runner local \
  --dry-run
```

再检查 `saturation`（需已有 `split_bams/per_spot_read_counts.tsv` 与 `coverage/host/**/*.CG.cov`）：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage saturation \
  --runner local \
  --dry-run
```

再检查 `summary`（需已有 `split_bams/per_spot_read_counts.tsv`、`demux/*.stats.json`、`qc/saturation/saturation_summary.tsv` 与 `coverage/` 下的 coverage 文件）：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage summary \
  --runner local \
  --dry-run
```

## 预期输出

`fastp_split` 后：

- `work/<sample>/shard_fastq/*.R1.fq.gz`
- `work/<sample>/shard_fastq/*.R2.fq.gz`
- `work/<sample>/shard_fastq/fastp.html`
- `work/<sample>/shard_fastq/fastp.json`

`demux_extract_bc` 后：

- `work/<sample>/demux/*.R1.demux.fq.gz`
- `work/<sample>/demux/*.R2.demux.fq.gz`
- `work/<sample>/demux/*.R1.spike-in.fq.gz`
- `work/<sample>/demux/*.R2.spike-in.fq.gz`
- `work/<sample>/demux/*.stats.json`

`align` 后：

- `work/<sample>/align_shards/*.cb.bam`
- `work/<sample>/align_shards/*.<spike_name>.bam`

`pool` 后：

- `work/<sample>/pooled/pooled.byCB.bam`
- 若配置了 `spike_in_index`：`work/<sample>/pooled/pooled.<spike_name>.sorted.bam`

`split` 后：

- `work/<sample>/split_bams/per_spot_read_counts.tsv`
- `work/<sample>/split_bams/**/*.sorted.bam`
- `work/<sample>/split_bams/**/*.sorted.bam.bai`

`mbias` 后（与 TAPS 相同路径；计数逻辑见 `scripts-emseq/mbias.py`）：

- `work/<sample>/qc/mbias/host.subsampled.sorted.bam`（及 `.bai`）
- `work/<sample>/qc/mbias/host.mbias.tsv`、`host.mbias.png`
- 若含 spike：`work/<sample>/qc/mbias/<spike_name>.mbias.tsv`、`.mbias.png`

`call` 后（兼容 TAPS coverage 契约）：

- `work/<sample>/pileup/**/*.vcf.gz`
- `work/<sample>/pileup/**/*.vcf.gz.tbi`
- `work/<sample>/coverage/host/<X_index>/<X_index>_<Y_index>.CG.cov`（已移除 `call_mito_chromosomes` 指定的 contig，默认不含 `chrM`）
- `work/<sample>/coverage/host_mito.CG.cov`（优先由 `qc/mbias/host.subsampled.sorted.bam` pileup 得到；否则由 per-spot coverage 汇总的线粒体位点得到，默认仅汇总 `chrM`）
- 若配置了 `spike_in_index`：`work/<sample>/coverage/<spike_name>.CG.cov`

`saturation` 后（与 TAPS 相同产物路径）：

- `work/<sample>/qc/saturation/saturation_curve.png`
- `work/<sample>/qc/saturation/saturation_summary.tsv`

`summary` 后（与 TAPS 相同产物路径）：

- `work/<sample>/summary/per_spot_summary.tsv`
- `work/<sample>/summary/sample_summary.tsv`
- `work/<sample>/summary/reads_heatmap.png`
- `work/<sample>/summary/cpg_site_count_heatmap.png`
- `work/<sample>/summary/mean_methylation_heatmap.png`
