# DBiT-Matrix

`DBiT-Matrix` 是面向 `DBiT-DNAme-TAPS` 的工作流：从原始 `FASTQ` 出发，完成条码提取、比对、spot 拆分、M-bias 质控、甲基化 calling，并输出 spot 级和 sample 级汇总结果。

- 当前版本：`1.2.0`
- 历史版本：`doc/log.md`

## 支持范围

- 当前协议：`DBiT-DNAme-TAPS`
- 固定主流程：`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> summary`
- 运行方式：`local`、`slurm`
- 环境管理：统一使用 `pixi`
- 配置入口：优先通过 `workflow/*.json` 管理
- 安全检查：所有 stage 支持 `--dry-run`

## 你会得到什么

最常用结果在 `work/<sample>/` 下：

- `summary/per_spot_summary.tsv`
- `summary/sample_summary.tsv`
- `summary/reads_heatmap.png`
- `summary/cpg_site_count_heatmap.png`
- `summary/mean_methylation_heatmap.png`
- `coverage/host/**/*.CG.cov`
- `coverage/host_mito.CG.cov`
- `coverage/<spike_name>.CG.cov`

`call` 阶段会按 batch 流式追加写出 `.CG.cov`，避免把整次 calling 结果长期累积在内存里。

## 你需要准备什么

- 双端测序数据：`R1 FASTQ`、`R2 FASTQ`
- barcode 白名单：`barcode1_whitelist`、`barcode2_whitelist`
- host 参考：`bwa_index`、`call_reference_file`
- 可选 spike-in 参考：`spike_in_index`
- workflow 配置文件（建议从 `workflow/dbit_taps_test.json` 复制）

## 三步快速开始

1) 安装环境

```bash
pixi install
```

2) 复制样例配置并修改最小必填项（`sample_id`、`r1`、`r2`、参考路径等）

3) 先 dry-run，再真实运行

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

建议按下面顺序检查：

1. `summary/sample_summary.tsv`（样本级总览）
2. `summary/per_spot_summary.tsv`（spot 级主表）
3. `summary/*.heatmap.png`（空间分布）

## 文档导航

首次用户建议阅读路径：
`README -> doc/setup.md -> doc/config.md -> doc/commands.md -> doc/outputs.md`

用户文档：

- `doc/setup.md`：环境和数据准备
- `doc/config.md`：`workflow/*.json` 配置说明（最小必改项）
- `doc/commands.md`：运行命令与 `local/slurm` 场景
- `doc/outputs.md`：主要结果文件与查看顺序
- `doc/stages.md`：各 stage 输入输出契约（参考手册）

维护与内部文档：

- `TEST.md`：回归/维护者检查
- `doc/progress.md`：内部里程碑与风险跟踪
- `doc/log.md`：版本变化记录
