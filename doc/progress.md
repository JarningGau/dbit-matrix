# Progress

## 当前状态

项目当前已闭环到 `mbias + call` MVP，实现了 host 抽样 M-bias QC、spike-in 全量 M-bias QC，以及 host spots 并行 calling 与 spike-in calling。

当前阶段链：

`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call`

## 已完成

- `fastp_split`
- `demux_extract_bc`
- `align`
- `pool`
- `split`
- `mbias`（MVP）
- `call`（MVP）

## 当前里程碑

`split -> mbias -> call` 端到端最小验收

目标：

- 增加 `split -> mbias -> call` 端到端验收命令
- 增加最小 smoke 级真实调用基线（输出文件数量与命名）

## 最近里程碑

- M1: `fastp_split + demux_extract_bc`
- M2: `align`
- M3: `pool`
- M4: `split`
- M5: `call` MVP
- M6: `mbias` MVP

## 风险点

- `call` 目前沿用旧 caller 的 pileup 逻辑，后续需评估性能与统计一致性
- `mbias` 当前仅输出 QC 表格，尚未自动反馈到 trimming/call 参数
- `host_mito` 当前是每 spot 一份输出，后续是否需要汇总层仍待讨论
- `split -> mbias -> call` 端到端验收尚未固化为稳定基线
