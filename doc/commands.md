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

- host: 从 `pooled/pooled.byCB.bam` 固定比例抽样（固定 seed）后 `sort + index + mbias`
- spike-in: 使用全量 `pooled/pooled.<spike_name>.sorted.bam` 直接做 mbias
- 输出到 `qc/mbias/`

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
- host 结果输出到 `coverage/host/` 与 `coverage/host_mito/`
- spike-in 结果输出到 `coverage/<spike_name>.CG.cov`

Slurm 下：

- 生成 `07_call_host.sbatch`（单作业内并行处理 spots）
- 生成 `07_call_spike_<spike_name>.sbatch`（每个 spike-in 一个 sbatch）
