
开发过程请遵循：
- 第一性原理
- 奥卡姆剃刀原则

## Rules（项目执行规范）

### 1) 迭代策略
- 先做最小可运行闭环（MVP），再做扩展，不做“先设计全量再实现”。
- 以阶段推进：`fastp_split -> demux_extract_bc -> align -> pool -> split -> call`。
- 每次迭代只引入一个主要变化点，避免多变量叠加导致排障困难。

### 2) 脚本与职责边界
- `scripts/*.py` 负责单步功能；`scripts/make_cmd.py` 负责命令生成与提交编排。
- 单步脚本保持薄封装，避免过度抽象；参数尽量显式、行为可预期。
- `--dry-run` 必须可用，作为所有变更的第一层验证入口。

### 3) 配置与命名
- 所有运行参数尽量配置化，优先放在 `workflow/*.json`。
- 配置文件命名必须语义化、稳定，避免频繁重命名造成文档和脚本漂移。
- stage 相关资源配置采用 step-specific 结构，不同阶段不要复用同一套 sbatch 参数。

### 4) HPC / Slurm 规范
- 不在 sbatch 模板中写与当前阶段无关的通用参数（例如不需要则不写 `--time`）。
- `demux_extract_bc` 在 Slurm 下采用“每个 chunk 一个 sbatch 脚本”，禁止单任务内 `for` 循环串行跑全量 chunk。
- `align` 在 Slurm 下同样采用“每个 chunk 一个 sbatch 脚本”，并在单个 chunk 脚本内按固定顺序执行该 chunk 的 spike-in 与 host 对齐。
- Slurm也用pixi来做环境管理，避免module load污染环境。

### 5) 输入输出契约
- 每个阶段必须定义稳定的输入/输出命名约定，并在 `README.md` 同步说明。
- 变更输出文件名或 reads header 格式时，必须同时更新文档与下游依赖。
- `demux` 输出需区分 matched 与 spike-in，统计文件需包含保留率和拒绝原因。
- `align` 输入约定：matched 使用 `*.R1/2.demux.fq.gz`，spike-in 使用 `*.R1/2.spike-in.fq.gz`。
- `align` 输出约定：host 输出 `<chunk>.cb.bam`；spike-in 输出 `<chunk>.<spike_name>.bam`。
- `align` 执行顺序约定：先跑 spike-in（可多个），再跑 host genome，确保行为可预测且便于回归。

### 5.1) Align Spike-in 配置约束
- `workflow/*.json` 中 align 相关参数保持显式：`bwa_index`、`bwa_threads`、`bwa_bin`、`sinto_bin`、`samtools_bin`。
- spike-in 参考配置统一使用 `spike_in_index`，支持对象或 `NAME=INDEX` 列表两种表示法，便于处理 0/1/N 个 spike-in 场景。

### 6) 环境与可复现性
- 统一使用 `pixi + uv` 管理环境，依赖以 `pixi.toml + pixi.lock` 为准。
- 文档与示例命令优先使用 `pixi run ...`，保证跨机器可复现。
- 新增依赖后必须更新 lockfile，并验证关键命令可在新锁定环境下执行。

### 7) 文档与测试
- 每次功能变更需同步更新 `README.md`（用法、输出、注意事项）与 `TEST.md`（最小验证步骤）。
- 测试的时候使用`pixi run ..`。
- 提交前至少执行：CLI `--help` 检查 + workflow dry-run 回归。

### 8) 提交与变更管理
- commit 信息描述“为什么做”，不是只写“改了什么”。
- 每次提交聚焦一个主题，避免将不相关改动混在同一 commit。
- 不回退或覆盖未明确授权的改动；遇到异常变更先确认再处理。

