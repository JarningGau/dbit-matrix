# DBiT-Matrix

空间组学数据分析工作流。当前支持 **TAPS**、**EMSeq** 与 **RNA**（MVP：`demux_extract_bc -> align`）入口，共享 `work/<sample>/` 目录结构与下游产物路径。

- 版本：`1.5.0`
- 变更记录：`doc/log.md`

## 入口选择

| 入口   | 编排脚本                | 配置示例                      |
|--------|------------------------|-------------------------------|
| TAPS   | `scripts/make_cmd.py`      | `workflow/dbit_taps_test.json`   |
| EMSeq  | `scripts-emseq/make_cmd.py` | `workflow/dbit_emseq_test.json`  |
| RNA    | `scripts-rna/make_cmd.py` | `workflow/dbit_rna_test.json` |

主流程：

`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> saturation -> summary`

`--stage all` 按序串联各 stage；所有 stage 支持 `--dry-run`。EMSeq / TAPS 另有实验性 `aggregate`（扁平汇总 per-spot `*.CG.cov`），需显式 `--stage aggregate`，不纳入 `--stage all`（见 `doc/stages.md` §11）。

## 主要产物

`work/<sample>/` 下的核心输出：

- `summary/sample_summary.tsv`、`summary/per_spot_summary.tsv`、`summary/*.heatmap.png`
- `qc/saturation/saturation_curve.png`、`qc/saturation/saturation_summary.tsv`
- `coverage/host/**/*.CG.cov`、`coverage/host_mito.CG.cov`、`coverage/<spike_name>.CG.cov`（可选）

字段说明与排错指南见 `doc/outputs.md`。

## 快速开始

### TAPS

1. `pixi install`
2. 复制 `workflow/dbit_taps_test.json`，配置必填字段：`sample_id`、`r1`、`r2`、参考路径、barcode 白名单
3. 验证配置后提交：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner local \
  --dry-run

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner local \
  --submit
```

Slurm 环境：`--runner local` 改为 `--runner slurm`。

### EMSeq

1. `pixi install`
2. 复制 `workflow/dbit_emseq_test.json`，按 `doc/emseq.md` 配置必填字段（含 `split_barcodes`、`call_reference_file`、`call_jobs`）
3. 验证配置后提交：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage all \
  --runner local \
  --dry-run

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage all \
  --runner local \
  --submit
```

Slurm 环境：`--runner local` 改为 `--runner slurm`；详见 `doc/emseq.md`。

### RNA (MVP)

1. `pixi install`
2. 复制 `workflow/dbit_rna_test.json` 并配置：`sample_id`、`r1`、`r2`、barcode 白名单、`linker_bc`、`umi_left`、`umi_len`、`star_genome_dir`、`gtf`
3. 先做 dry-run 再提交：

```bash
pixi run python scripts-rna/make_cmd.py \
  --workflow-config workflow/dbit_rna_test.json \
  --stage all \
  --runner local \
  --dry-run

pixi run python scripts-rna/make_cmd.py \
  --workflow-config workflow/dbit_rna_test.json \
  --stage all \
  --runner local \
  --submit
```

Slurm 环境：`--runner local` 改为 `--runner slurm`。

RNA MVP 说明：`align` 使用 STARsolo（10xv2-like `CB_UMI_Simple`），输出为 `work/<sample>/solo/` 下的 matrix 结果，并在 `Gene/raw` 与 `Gene/filtered` 目录写出 `barcodes_pos.tsv`（列 `barcode\\tx\\ty`，0-based）；RNA 入口不使用 chunk 概念。

## 结果检查顺序

1. `summary/sample_summary.tsv`
2. `summary/per_spot_summary.tsv`
3. `summary/*.heatmap.png`
4. `qc/saturation/saturation_summary.tsv`

## 文档索引

### TAPS

`README` → `doc/setup.md` → `doc/config.md` → `doc/commands.md` → `doc/outputs.md` → `doc/stages.md`

- `doc/setup.md`：环境配置与测试数据
- `doc/config.md`：配置文件必填字段
- `doc/commands.md`：运行命令（`local` / `slurm`）
- `doc/outputs.md`：结果解读
- `doc/stages.md`：各 stage 输入输出契约

### EMSeq

- `doc/emseq.md`：EMSeq 配置与使用

### 维护者

- `TEST.md`：TAPS 回归测试
- `TEST-emseq.md`：EMSeq 回归测试
- `TEST-rna.md`：RNA 回归测试

### 内部

- `doc/progress.md`：里程碑与风险
- `doc/log.md`：版本记录
- `doc/TODO.md`：待办事项
