# Commands

本页面向 1.1.2 用户，只保留 `scripts/make_cmd.py` 的常用命令。默认推荐先 `--dry-run`，确认输入路径、stage 展开顺序和 runner 设置都正确。

## 最常用

查看帮助：

```bash
pixi run python scripts/make_cmd.py --help
```

完整流程 dry-run：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner local \
  --dry-run
```

完整流程本地执行：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner local \
  --submit
```

完整流程提交到 Slurm：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner slurm \
  --submit
```

说明：

- 推荐显式传 `--stage all`
- `all` 的固定展开顺序是 `fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> summary`
- 未显式传 `--stage` 时，会优先读取 `workflow/*.json` 里的 `stage`

## Runner 说明

`local`：

- 生成每个 stage 的 `.sh`
- 额外生成一个串行入口 `run.sh`
- 使用 `--submit` 时执行 `run.sh`

`slurm`：

- 生成每个 stage 的 `.sbatch`
- 额外生成一个入口 `run.sbatch`
- 使用 `--submit` 时提交 `run.sbatch`
- stage 依赖通过 `afterok` 管理

## 按 Stage 运行

如果你只想检查某一段流程，可以显式指定 `--stage`。

### `demux_extract_bc`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage demux_extract_bc \
  --runner local \
  --dry-run
```

用途：

- 从 `shard_fastq/` 读取 chunk FASTQ
- 提取 barcode
- 输出 matched reads、spike-in reads 和统计信息

### `align`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage align \
  --runner local \
  --dry-run
```

用途：

- 读取 `demux/` 输出
- 先对 spike-in 比对，再对 host 比对
- 生成 host 和 spike-in 的 BAM 分片

### `pool`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage pool \
  --runner local \
  --dry-run
```

用途：

- 合并 host BAM 分片
- 合并并排序 spike-in BAM 分片

### `split`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage split \
  --runner local \
  --dry-run
```

用途：

- 读取 `pooled/pooled.byCB.bam`
- 按 `CB:Z:<x>+<y>` 拆成 spot BAM
- 对 spot BAM 做后续排序

补充：

- `split_smoke=true` 或 `--split-smoke` 时，只输出最多 16 个非空 spot

### `mbias`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage mbias \
  --runner local \
  --dry-run
```

用途：

- 对 host 和/或 spike-in 生成 M-bias 质控结果
- 输出 `qc/mbias/*.mbias.tsv` 和 `qc/mbias/*.mbias.png`
- host 抽样 BAM 会固定输出到 `qc/mbias/host.subsampled.sorted.bam`，供 `call` 复用

### `call`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage call \
  --runner local \
  --call-mode all \
  --dry-run
```

用途：

- 对 host spot BAM 进行 calling
- 生成聚合 `host_mito` calling 结果
- 对 spike-in 生成聚合 calling 结果

补充：

- `--call-mode all`：同时生成 host 和 spike-in
- `--call-mode host`：只生成 host 结果
- `--call-mode spike`：只生成 spike-in 结果

### `summary`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage summary \
  --runner local \
  --dry-run
```

用途：

- 读取 `coverage/` 和 `split_bams/per_spot_read_counts.tsv`
- 输出 `summary/per_spot_summary.tsv`
- 输出 `summary/sample_summary.tsv`
- 输出 `summary/reads_heatmap.png`
- 输出 `summary/cpg_site_count_heatmap.png`
- 输出 `summary/mean_methylation_heatmap.png`

## 常见组合

只验证 Slurm 命令是否能正确展开：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner slurm \
  --dry-run
```

只检查 `mbias`：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage mbias \
  --runner local \
  --dry-run
```

只检查 `call` 的 host 结果：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage call \
  --runner local \
  --call-mode host \
  --dry-run
```
