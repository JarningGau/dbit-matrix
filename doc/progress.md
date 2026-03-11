# Progress

## 当前状态

项目当前已闭环到 `split`，尚未进入 `call` 的实现闭环。

当前阶段链：

`fastp_split -> demux_extract_bc -> align -> pool -> split -> call`

## 已完成

- `fastp_split`
- `demux_extract_bc`
- `align`
- `pool`
- `split`

## 当前里程碑

`call`

目标：

- 明确 `call` 的输入输出契约
- 增加单步脚本与 `make_cmd.py --stage call`
- 补齐本地 / Slurm dry-run
- 补齐 `README.md` 与 `TEST.md`

## 最近里程碑

- M1: `fastp_split + demux_extract_bc`
- M2: `align`
- M3: `pool`
- M4: `split`

## 风险点

- `call` 尚未确定最终输出格式
- 当前用户文档主要围绕 `DBiT-DNAme-TAPS`
- 端到端验收目前停在 `split`
