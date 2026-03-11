# Progress

## 当前状态

项目当前已闭环到 `call` MVP，实现了 host spots 并行 calling 与 spike-in calling。

当前阶段链：

`fastp_split -> demux_extract_bc -> align -> pool -> split -> call`

## 已完成

- `fastp_split`
- `demux_extract_bc`
- `align`
- `pool`
- `split`
- `call`（MVP）

## 当前里程碑

`split -> call` 端到端最小验收

目标：

- 增加 `split -> call` 端到端验收命令
- 增加最小 smoke 级真实调用基线（输出文件数量与命名）

## 最近里程碑

- M1: `fastp_split + demux_extract_bc`
- M2: `align`
- M3: `pool`
- M4: `split`
- M5: `call` MVP

## 风险点

- `call` 目前沿用旧 caller 的 pileup 逻辑，后续需评估性能与统计一致性
- `host_mito` 当前是每 spot 一份输出，后续是否需要汇总层仍待讨论
- `split -> call` 端到端验收尚未固化为稳定基线
