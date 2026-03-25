# TEST-RNA (Maintainer)

DBiT-RNA 独立入口维护者回归。输入输出契约见 `doc/stages.md`（RNA 小节）。

## 范围

- 入口：`scripts-rna/make_cmd.py`
- Stages：`demux_extract_bc`、`align`，以及 `all`（生成 `commands/run.sh` / `run.sbatch` 串联两阶段）
- 单步脚本：`scripts-rna/extract_bc.py`、`scripts-rna/align.py`
- 样例配置：`workflow/dbit_rna_test.json`

## 最小验收

提交前至少完成：

1. CLI `--help` 可用
2. 两个具名 stage 在 `local` 与 `slurm` 下均能 dry-run；`--stage all` 在 `local` 与 `slurm` 下 dry-run
3. 用样例配置完成一次本地真实生成（可逐 stage，也可 `--stage all` 生成 `run.sh` 后执行）

### CLI

```bash
pixi run python scripts-rna/make_cmd.py --help
pixi run python scripts-rna/extract_bc.py --help
pixi run python scripts-rna/align.py --help
```

### Dry-run：`local`（两阶段 + `all`）

```bash
CFG=workflow/dbit_rna_test.json
STAGES="demux_extract_bc align"
for stage in $STAGES; do
  pixi run python scripts-rna/make_cmd.py \
    --workflow-config "$CFG" \
    --stage "$stage" \
    --runner local \
    --dry-run
done

pixi run python scripts-rna/make_cmd.py \
  --workflow-config "$CFG" \
  --stage all \
  --runner local \
  --dry-run
```

### Dry-run：`slurm`（两阶段 + `all`）

将上述命令中 `--runner local` 改为 `--runner slurm`。

### 通过标准（dry-run）

- 命令正常展开，无缺失参数、无路径解析错误
- `demux_extract_bc` 调用 `scripts-rna/extract_bc.py`
- `align` 调用 `scripts-rna/align.py`
- 输出目录指向 `demux`、`solo`、`commands` 等预期位置
- `stage=all` 在 dry-run 下不落盘
- `runner=slurm` 且**非** `--dry-run` 生成脚本时：各 `.sbatch` 含 `#SBATCH --output` / `--error`，默认 `work/<sample>/logs/<stem>_%x_%j.out` / `.err`（`demux_extract_bc`、`align`；`run.sbatch` driver 为 `rna_run_`）；需要时校验 `logs/` 已创建

## 真实生成 Smoke（本地）

```bash
rm -rf work/test-DBiT-RNA

CFG=workflow/dbit_rna_test.json
STAGES="demux_extract_bc align"
for stage in $STAGES; do
  pixi run python scripts-rna/make_cmd.py \
    --workflow-config "$CFG" \
    --stage "$stage" \
    --runner local \
    --submit
done
```

或一次性串联执行：

```bash
rm -rf work/test-DBiT-RNA
pixi run python scripts-rna/make_cmd.py \
  --workflow-config workflow/dbit_rna_test.json \
  --stage all \
  --runner local \
  --submit
```

### 输出检查点

Smoke 通过时至少存在（路径相对于 `work/<sample>/`）：

- `demux/<sample>.R1.clean.fq.gz`、`demux/<sample>.R2.clean.fq.gz`、`demux/<sample>.stats.json`
- `solo/`（包含 STARsolo matrix 产物，如 `solo/star.Solo.out/Gene/`）
- `solo/star.Solo.out/Gene/raw/barcodes_pos.tsv`、`solo/star.Solo.out/Gene/filtered/barcodes_pos.tsv`
- `commands/01_demux_extract_bc.sh`、`commands/02_align.sh`（按 `--stage all` 运行时）
- `commands/run.sh`（按 `--stage all` 运行时）

### 关键行为检查

- `demux_extract_bc`：R1 结构按 `BC2-linker_bc-BC1-UMI_left-UMI` 解析
- `demux_extract_bc`：R1 clean 读长应为 `8+8+umi_len`（样例配置为 26）
- `demux_extract_bc`：read name 包含 `BC2+BC1:UMI:original_name`
- `align`：输入为 `demux/<sample>.R1.clean.fq.gz` 与 `demux/<sample>.R2.clean.fq.gz`（无 chunk）
- `align`：运行 STARsolo 并将 matrix 结果写入 `solo/`
- `align`：生成 `Gene/raw` 与 `Gene/filtered` 的 `barcodes_pos.tsv`（列：`barcode\tx\ty`，`x/y` 为 whitelist 行号 0-based）

## 数据说明

回归前确认 `workflow/dbit_rna_test.json` 中 `r1` / `r2`、`star_genome_dir`、`gtf` 路径存在；并确认 barcode 白名单、`linker_bc`、`umi_left`、`umi_len` 与实验设计一致。
