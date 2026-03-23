
开发过程请遵循：
- 第一性原理
- 奥卡姆剃刀原则

## 当前阶段

- 主线固定（EMSeq 独立入口）：`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> saturation -> summary`
- 当前已闭环（EMSeq MVP）：`fastp_split`、`demux_extract_bc`、`align`、`pool`、`split`、`mbias`、`call`、`saturation`、`summary`
- EMSeq 编排入口支持 `--stage all`：生成 `work/<sample>/commands/run.sh` 或 `run.sbatch`，按主线顺序串联各 stage（`--dry-run` 不落盘）。
- 当前主里程碑：`split -> mbias -> call -> saturation -> summary` 端到端最小闭环验收
- 说明：`mbias` 单步实现为 `scripts-emseq/mbias.py`（与 TAPS 的 `scripts/mbias.py` 化学判定不同）。

## 必须遵守

- 先做 MVP，再做扩展；不要先设计全量方案。
- 每次迭代只引入一个主要变化点。
- `scripts/*.py` 与 `scripts-emseq/*.py` 负责单步功能；workflow 编排入口：TAPS 用 `scripts/make_cmd.py`，EMSeq 用 `scripts-emseq/make_cmd.py`。
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
- `mbias`：按 host/spike 拆分 sbatch；资源配置使用 `slurm.mbias.host` / `slurm.mbias.spike`。
- `call`：按 host/spike 拆分 sbatch（host 1 个；spike 按 `spike_name` 可为多个），资源配置保持 step-specific（如 `slurm.call.host` / `slurm.call.spike`）。
- `saturation`：单 job 一个 sbatch；资源配置使用 `slurm.saturation`。
- `summary`：单 job 一个 sbatch；资源配置使用 `slurm.summary`。
- stage 资源配置采用 step-specific 结构，不同阶段不要混用一套参数。

## 关键契约

- 文档里的输入输出契约必须与代码一致；改契约时同步更新 `README.md` 与下游依赖。
- `demux`：输出 host demux 与 spike-in FASTQ（`*.demux.fq.gz` 与 `*.spike-in.fq.gz`）；统计文件必须包含保留率和拒绝原因。
- `align`：host 输入 `*.demux.fq.gz`，spike-in 输入 `*.spike-in.fq.gz`；输出分别为 `<chunk>.cb.bam` 与 `<chunk>.<spike_name>.bam`。
- `align`：`spike_in_index` 统一支持对象或 `NAME=INDEX` 列表。
- `pool`：host 最终输出 `pooled/pooled.byCB.bam`；spike-in 输出 `pooled/pooled.<spike_name>.sorted.bam`。
- `split`：输入 `pooled/pooled.byCB.bam`，按 `CB:Z:<x>+<y>` 解析 spot。
- `split`：smoke 模式最多输出 16 个非空 spot BAM。

## 变更前后检查

- 功能变更后同步更新 `README.md` 与对应回归文档（TAPS：`TEST.md`；EMSeq：`TEST-emseq.md`）。
- 推进一个阶段后，同步更新 README 中的“当前进度 / 当前里程碑”。
- 提交前至少执行：CLI `--help` 检查 + workflow dry-run 回归。

## 提交规则

- commit 信息写“为什么做”，不要只写“改了什么”。
- 每次提交只聚焦一个主题。
- 建议主题行（subject）格式统一为“范围/阶段 + 为什么（一句话）”，优先用当前工作流阶段做前缀，便于回溯与检索：
  - `EMSeq <stage>: <why>`（例如 `EMSeq call: ...` / `EMSeq split: ...`）
  - 或 `docs: <why>`、`workflow(emseq): <why>`、`scripts-emseq: <why>` 等（范围清晰即可）
- subject 避免空泛词（如 update/refine/tweak）；尽量写清“要解决的问题/要保障的约束”，例如：
  - “preserve downstream BAM contracts / keep dry-run shell-safe / keep resource settings step-specific”
- 当同一提交同时涉及代码 + workflow 配置 + 文档时，subject 只写**主因**；其余放 body，用于说明“哪些契约被保持/同步了”。
- body（可选但推荐）用最少文字回答三件事：
  - 为什么要改（动机/问题）
  - 不变量/约束（哪些输入输出契约不能破；性能/资源边界）
  - 如何验证（至少：CLI `--help` + workflow `--dry-run`）
- 任何“输入输出契约变化”必须在 body 显式写出，并同步更新 `README.md` 与下游依赖；若为破坏性变更，明确标注 `BREAKING:`。
- 避免噪音行淹没关键信息（例如与变更动机无关的自动署名行）；commit 信息以“可复用的决策记录”为准。
- 不回退或覆盖未明确授权的改动；遇到异常先确认。

