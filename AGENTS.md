
开发过程请遵循：
- 第一性原理
- 奥卡姆剃刀原则

## 当前阶段

- 主线固定：`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call`
- 当前已闭环：`fastp_split`、`demux_extract_bc`、`align`、`pool`、`split`、`mbias`、`call`（MVP）
- 当前主里程碑：`split -> mbias -> call` 端到端最小闭环验收

## 必须遵守

- 先做 MVP，再做扩展；不要先设计全量方案。
- 每次迭代只引入一个主要变化点。
- `scripts/*.py` 负责单步功能；`scripts/make_cmd.py` 负责 workflow 编排。
- 单步脚本保持薄封装，参数显式、行为可预测。
- 所有阶段必须支持 `--dry-run`。
- 所有运行参数优先放到 `workflow/*.json`。
- 新增依赖后必须更新 `pixi.lock`。

## Slurm 规则

- Slurm 使用 `pixi` 环境，不依赖 `module load`。
- 不在 sbatch 模板中写与当前阶段无关的通用参数。
- `demux_extract_bc`：每个 chunk 一个 sbatch。
- `align`：每个 chunk 一个 sbatch；单脚本内顺序执行 spike-in -> host。
- `pool`：拆成两个 job，一个 spike-in，一个 host。
- `split`：`split_bams` 与 `sort` 分成两个 sbatch，并用依赖关系串联。
- stage 资源配置采用 step-specific 结构，不同阶段不要混用一套参数。

## 关键契约

- 文档里的输入输出契约必须与代码一致；改契约时同步更新 `README.md` 与下游依赖。
- `demux`：输出 matched 与 spike-in；统计文件必须包含保留率和拒绝原因。
- `align`：host 输入 `*.demux.fq.gz`，spike-in 输入 `*.spike-in.fq.gz`；输出分别为 `<chunk>.cb.bam` 与 `<chunk>.<spike_name>.bam`。
- `align`：`spike_in_index` 统一支持对象或 `NAME=INDEX` 列表。
- `pool`：host 最终输出 `pooled/pooled.byCB.bam`；spike-in 输出 `pooled/pooled.<spike_name>.sorted.bam`。
- `split`：输入 `pooled/pooled.byCB.bam`，按 `CB:Z:<x>+<y>` 解析 spot。
- `split`：smoke 模式最多输出 16 个非空 spot BAM。

## 变更前后检查

- 功能变更后同步更新 `README.md` 与 `TEST.md`。
- 推进一个阶段后，同步更新 README 中的“当前进度 / 当前里程碑”。
- 提交前至少执行：CLI `--help` 检查 + workflow dry-run 回归。

## 提交规则

- commit 信息写“为什么做”，不要只写“改了什么”。
- 每次提交只聚焦一个主题。
- 不回退或覆盖未明确授权的改动；遇到异常先确认。

