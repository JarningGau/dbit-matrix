# Progress

内部里程碑与风险跟踪。对外能力说明以 `README.md` 为准。

## 当前状态

EMSeq 独立入口主线已闭环：

`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> saturation -> summary`

DBiT-RNA MVP 已接入独立入口（当前范围）：

`demux_extract_bc -> align`

RNA 当前 `align` 已从占位 BAM 迁移为 STARsolo matrix-only 输出（`work/<sample>/solo/`），并采用非 chunk 模式。

## 已完成（主线）

- `fastp_split`
- `demux_extract_bc`
- `align`
- `pool`
- `split`
- `mbias`
- `call`
- `saturation`
- `summary`

## 后续工作（非阻塞）

- `call`：pileup 性能与统计一致性持续评估
- `mbias`：当前仅输出 QC，是否自动反馈 trimming/call 参数待产品决策
- `host_mito`：单文件聚合是否扩展为分层统计待讨论
- 端到端 smoke 基线与资源经验值持续沉淀（见 `TEST-emseq.md`）

## 历史里程碑

- M1: `fastp_split + demux_extract_bc`
- M2: `align`
- M3: `pool`
- M4: `split`
- M5–M7: `call` / `mbias` / `summary` MVP
- 补齐 `saturation` 与全链路汇总

## 风险点

与「后续工作」一致，不重复展开。
