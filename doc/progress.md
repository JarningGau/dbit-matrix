# Progress

本页用于内部里程碑和风险跟踪，不作为用户支持范围说明。  
对外能力与使用入口以 `README.md` 为准。

## 当前状态

项目当前已闭环到 `summary`，在既有 `mbias + call` MVP 基础上补齐了结果汇总与可视化出口。

当前阶段链：

`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> summary`

## 已完成

- `fastp_split`
- `demux_extract_bc`
- `align`
- `pool`
- `split`
- `mbias`（MVP）
- `call`（MVP）

## 当前里程碑

`split -> mbias -> call -> summary` 端到端最小验收

目标：

- 增加 `split -> mbias -> call -> summary` 端到端验收命令
- 增加最小 smoke 级真实调用基线（输出文件数量与命名）

## 最近里程碑

- M1: `fastp_split + demux_extract_bc`
- M2: `align`
- M3: `pool`
- M4: `split`
- M5: `call` MVP
- M6: `mbias` MVP
- M7: `summary` MVP

## 风险点

- `call` 目前沿用旧 caller 的 pileup 逻辑，后续需评估性能与统计一致性
- `mbias` 当前仅输出 QC 表格，尚未自动反馈到 trimming/call 参数
- `host_mito` 当前为单个聚合输出，后续是否需要增加分层统计仍待讨论
- `split -> mbias -> call -> summary` 端到端验收尚未固化为稳定基线
