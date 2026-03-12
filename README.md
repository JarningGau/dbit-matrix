# DBiT-Matrix

`DBiT-Matrix` 是一套面向 DBiT-Omics 上游处理的 workflow，目标是同时支持本地运行与 HPC/Slurm 提交，并保持环境可复现。

当前 README 只保留总览；详细说明按主题拆到 `doc/` 下，采用渐进披露方式组织。

## 当前范围

- 当前重点记录：`DBiT-DNAme-TAPS`
- 已闭环阶段：`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call`
- 下一主里程碑：`split -> mbias -> call` 端到端最小验收

## 快速开始

安装环境：

```bash
pixi lock
pixi install
```

最小 dry-run：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --dry-run
```

## 当前进度

| 阶段 | 状态 | 说明 |
| --- | --- | --- |
| `fastp_split` | 已完成 | 支持本地与 Slurm |
| `demux_extract_bc` | 已完成 | 支持 matched/spike-in 分流与统计 |
| `align` | 已完成 | 支持 host 与多 spike-in |
| `pool` | 已完成 | 支持 host/spike-in 汇总与排序 |
| `split` | 已完成 | 支持按 spot 拆分、smoke、后续并行排序 |
| `mbias` | 已完成（MVP） | 支持 host 抽样 + sort/index 与 spike-in 全量、基于真实 CpG 位点的 M-bias QC；host 抽样 BAM 可被 `call` 复用 |
| `call` | 已完成（MVP） | 支持 `call_mode=all/host/spike`、host spots 并行 calling、spike-in calling，以及聚合 `host_mito` 调用（优先复用 `qc/mbias/host.subsampled.sorted.bam`） |

## 里程碑

- M1: `fastp_split + demux_extract_bc`
- M2: `align`
- M3: `pool`
- M4: `split`
- M5: `call` MVP（已完成）
- M6: `mbias` MVP（已完成，位于 `call` 前）
- M7: `split -> mbias -> call` 端到端最小闭环验收

## 文档导航

- `doc/setup.md`: 环境、测试数据、样例配置
- `doc/stages.md`: 各阶段输入、输出与行为约定
- `doc/commands.md`: `make_cmd.py` 的本地/Slurm 命令生成方式
- `doc/progress.md`: 当前开发进度、已完成里程碑、下一步
- `TEST.md`: 最小回归步骤
- `AGENTS.md`: 开发规则与阶段约束

## 目录概览

- `scripts/`: 单步脚本与 workflow 编排入口
- `workflow/`: 参数与 Slurm 资源配置
- `configs/`: whitelist 等静态配置
- `data/raw/`: 原始数据
- `work/`: 工作目录与阶段产物

## 备注

- 环境统一使用 `pixi`
- 如需引入 PyPI 包，优先使用 `pixi add --pypi ...`
- 文档示例命令统一使用 `pixi run ...`
