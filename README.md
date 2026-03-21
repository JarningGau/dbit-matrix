# DBiT-Matrix

`DBiT-Matrix` 是面向空间 DNA 甲基化数据的工作流集合。当前正式主线是 `DBiT-DNAme-TAPS`；`EMSeq` 作为独立入口保留在试验阶段。

- 当前版本：`1.4.0`
- 版本变化：`doc/log.md`

## 当前支持范围

### TAPS 主线

- 协议：`DBiT-DNAme-TAPS`
- 固定主流程：`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> saturation -> summary`
- 入口脚本：`scripts/make_cmd.py`
- 运行方式：`local`、`slurm`
- 环境管理：统一使用 `pixi`
- 配置入口：优先通过 `workflow/*.json`
- 安全检查：所有 stage 支持 `--dry-run`

### EMSeq 独立入口

- 当前状态：试验性 / MVP
- 当前支持：`fastp_split -> demux_extract_bc -> align -> pool -> split -> call`
- 入口脚本：`scripts-emseq/make_cmd.py`
- 说明页面：`doc/emseq.md`
- `call` 阶段默认将 `chrM` 从 `coverage/host/**/*.CG.cov` 中移除，并汇总为 `coverage/host_mito.CG.cov`

## 你会得到什么

TAPS 主线最常用结果位于 `work/<sample>/`：

- `summary/sample_summary.tsv`
- `summary/per_spot_summary.tsv`
- `summary/reads_heatmap.png`
- `summary/cpg_site_count_heatmap.png`
- `summary/mean_methylation_heatmap.png`
- `qc/saturation/saturation_curve.png`
- `qc/saturation/saturation_summary.tsv`
- `coverage/host/**/*.CG.cov`
- `coverage/host_mito.CG.cov`
- `coverage/<spike_name>.CG.cov`

`summary/sample_summary.tsv` 当前除甲基化汇总外，还包含 `raw_reads`、`barcoded_reads`、`host_mapped_reads`、`host_valid_reads`、`<spike_name>_mapped_reads`、`barcoded_reads_rate`、`valid_reads_rate` 和 `saturation_rate`。

## 三步快速开始

1. 安装环境：

```bash
pixi install
```

2. 复制 `workflow/dbit_taps_test.json`，先只修改最小必填字段：`sample_id`、`r1`、`r2`、参考路径、barcode 白名单。

3. 先 dry-run，再真实运行：

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

如需提交到 Slurm，把 `--runner local` 改为 `--runner slurm`。

## 跑完后先看什么

建议先看：

1. `summary/sample_summary.tsv`
2. `summary/per_spot_summary.tsv`
3. `summary/*.heatmap.png`
4. `qc/saturation/saturation_summary.tsv`

## 文档导航

首次使用 TAPS，建议阅读：
`README -> doc/setup.md -> doc/config.md -> doc/commands.md -> doc/outputs.md -> doc/stages.md`

### TAPS 用户

- `doc/setup.md`：环境、测试数据和准备材料
- `doc/config.md`：`workflow/*.json` 最小必改项与常用字段
- `doc/commands.md`：运行命令与 `local/slurm` 场景
- `doc/outputs.md`：主要产物与查看顺序
- `doc/stages.md`：各 stage 输入输出契约

### EMSeq 用户

- `doc/emseq.md`：EMSeq 当前支持范围、最小配置和首次运行命令

### 维护者 / 开发者

- `TEST.md`：TAPS 回归检查
- `TEST-emseq.md`：EMSeq 回归检查

### 内部材料

- `doc/progress.md`：里程碑、风险和下一步
- `doc/log.md`：版本变化记录
- `doc/TODO.md`：设计备忘和待决策事项
