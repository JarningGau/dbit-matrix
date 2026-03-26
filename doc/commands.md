# Run Commands

本页是 `scripts/make_cmd.py` 的运行手册，只覆盖用户侧运行命令与 `runner` 行为。回归型检查命令请看 `TEST.md`；EMSeq 独立入口请看 `doc/emseq.md`。

## 最常用命令

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

完整流程 Slurm 执行：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner slurm \
  --submit
```

## 固定行为

- 推荐显式传 `--stage all`
- `all` 固定展开：`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> saturation -> summary`
- 未显式传 `--stage` 时，会优先读取 `workflow/*.json` 中的 `stage`
- 全部 stage 都支持 `--dry-run`

## Runner 说明

`local`：

- 生成每个 stage 的 `.sh`
- 额外生成串行入口 `run.sh`
- `--submit` 时执行 `run.sh`

`slurm`：

- 生成每个 stage 或 chunk 的 `.sbatch`
- 额外生成入口 `run.sbatch`
- `--submit` 时由客户端提交依赖 DAG
- stage 间依赖统一用 `afterok`
- 不依赖计算节点内 nested `sbatch`

`split` 在 `slurm` 下固定拆成两段：

- `05_split_bams.sbatch`
- `05_split_sort.sbatch`

并通过 `05_split_submit.sh` 串联为
`05_split_bams.sbatch -> afterok -> 05_split_sort.sbatch`。

## 按 Stage dry-run

### `demux_extract_bc`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage demux_extract_bc \
  --runner local \
  --dry-run
```

### `align`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage align \
  --runner local \
  --dry-run
```

### `pool`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage pool \
  --runner local \
  --dry-run
```

### `split`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage split \
  --runner local \
  --dry-run
```

### `mbias`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage mbias \
  --runner local \
  --dry-run
```

### `call`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage call \
  --runner local \
  --call-mode all \
  --dry-run
```

### `saturation`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage saturation \
  --runner local \
  --dry-run
```

### `summary`

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage summary \
  --runner local \
  --dry-run
```

### `aggregate`（实验性）

不纳入 `all`；需显式指定 stage。

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage aggregate \
  --runner local \
  --dry-run
```

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage aggregate \
  --runner slurm \
  --dry-run
```

### `methscan_prepare`（实验性）

不纳入 `all`；需显式指定 stage。依赖已存在的 `coverage/host/**/*.CG.cov`。

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage methscan_prepare \
  --runner local \
  --dry-run
```

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage methscan_prepare \
  --runner slurm \
  --dry-run
```

## 常见组合

只验证 Slurm 命令是否能正确展开：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner slurm \
  --dry-run
```

只检查 host calling：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage call \
  --runner local \
  --call-mode host \
  --dry-run
```

只检查 spike calling：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage call \
  --runner local \
  --call-mode spike \
  --dry-run
```
