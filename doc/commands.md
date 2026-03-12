# Commands

本页只记录 workflow 入口 `scripts/make_cmd.py` 的常用命令。

## 通用

最小 dry-run：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --dry-run
```

直接执行：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --submit
```

切换到 Slurm：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --runner slurm \
  --submit
```

提示：

- 不显式传 `--stage` 时，按 `workflow/*.json` 中的 `stage` 执行（若缺省则是 `fastp_split`）
- 需要 workflow 一键入口时，请显式传 `--stage all`

## `all`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner local \
  --dry-run
```

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner slurm \
  --dry-run
```

说明：

- 展开顺序固定：`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> summary`
- local 下会生成每个 stage 的 `.sh`，并额外生成 `run.sh`
- local + `--submit`：提交 `run.sh`
- slurm 下会生成每个 stage 的 `.sbatch`，并额外生成 `run.sbatch`
- slurm + `--submit`：提交 `run.sbatch`；由 `run.sbatch` 调用 `sbatch --dependency=afterok:...` 管理阶段依赖

## `demux_extract_bc`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage demux_extract_bc \
  --dry-run
```

说明：

- 扫描 `shard_fastq/*.R1.fq.gz`
- Slurm 下每个 chunk 一个 sbatch

## `align`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage align \
  --dry-run
```

说明：

- 扫描 `demux/*.R1.demux.fq.gz`
- 本地模式按 chunk 顺序执行
- Slurm 下每个 chunk 一个 sbatch

## `pool`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage pool \
  --dry-run
```

说明：

- 扫描 `align_shards/*.cb.bam` 与 `align_shards/*.<spike_name>.bam`
- 本地模式先 spike-in 后 host
- Slurm 下生成 `04_pool_spike.sbatch` 与 `04_pool_host.sbatch`

## `split`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage split \
  --dry-run
```

说明：

- 读取 `pooled/pooled.byCB.bam`
- 先调用 `scripts/split_bams.py`
- 再调用 `scripts/bam_sort_parallel.py`
- `split_smoke=true` 或 `--split-smoke` 时仅输出最多 16 个非空 spot

Slurm 下：

- 生成 `05_split_bams.sbatch`
- 生成 `05_split_sort.sbatch`
- `sort` 通过 `afterok` 依赖 `split_bams`

## `mbias`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage mbias \
  --dry-run
```

说明：

- `mbias` 默认 `mode=spike`（只分析 spike-in）
- host: 从 `pooled/pooled.byCB.bam` 固定比例抽样（内部固定 seed）后 `sort + index + mbias`
- host: 抽样 BAM 固定输出到 `qc/mbias/host.subsampled.sorted.bam`，可被 `call` 聚合 `host_mito` 复用
- host: `R1` 按右对齐 cycle 统计（适配 demux 的 left trimming）
- spike-in: 使用全量 `pooled/pooled.<spike_name>.sorted.bam` 直接做 mbias
- host 使用 `call_reference_file`，spike-in 使用 `spike_in_index`
- 甲基化率只在参考序列真实 `CpG` 位点上统计：`(TG+CA)/(TG+CA+CG)`
- 输出到 `qc/mbias/`，每个样本同时产出 `*.mbias.tsv` 与 `*.mbias.png`

Slurm 下：

- 生成 `06_mbias_host.sbatch`
- 生成 `06_mbias_spike_<spike_name>.sbatch`

## `call`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage call \
  --dry-run
```

说明：

- 扫描 `split_bams/**/*.sorted.bam` 与 `pooled/pooled.<spike_name>.sorted.bam`
- 本地模式调用 `scripts/call.py`，在单脚本中并行 host spots
- host spot 结果输出到 `coverage/host/`
- `host_mito` 输出为单个聚合文件 `coverage/host_mito.CG.cov`
- `host_mito` 优先复用 `qc/mbias/host.subsampled.sorted.bam`，缺失时自动从 `pooled/pooled.byCB.bam` 抽样并排序
- spike-in 结果输出到 `coverage/<spike_name>.CG.cov`

Slurm 下：

- 生成 `07_call_host.sbatch`（单作业内并行处理 spots）
- 生成 `07_call_spike_<spike_name>.sbatch`（每个 spike-in 一个 sbatch）

## `summary`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage summary \
  --dry-run
```

说明：

- 读取 `coverage/` 与 `split_bams/per_spot_read_counts.tsv`
- 输出 `summary/per_spot_summary.tsv` 与 `summary/sample_summary.tsv`
- 缺失输入保持固定 summary 列并写 `NA`

Slurm 下：

- 生成 `08_summary.sbatch`
