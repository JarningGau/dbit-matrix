# Version Log

版本变更记录。当前支持范围与入口以 `README.md` 为准。

## 1.5.0

- 新增 EMSeq 独立入口：`scripts-emseq/make_cmd.py`，主线闭环为 `fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> saturation -> summary`
- 新增 EMSeq 单步脚本：`extract_bc.py`、`aligner.py`、`mbias.py`、`call.py`，并补齐 `workflow/dbit_emseq_test.json` / `workflow/dbit_emseq_slurm.json`
- `--stage all` 支持在 `local` / `slurm` 下按主线顺序串联全流程；继续保持各 stage `--dry-run` 与下游产物路径契约
- 新增 EMSeq 用户文档与维护者回归：`doc/emseq.md`、`TEST-emseq.md`，并同步 README 与 stage/config/output 文档

## 1.4.0

- 文档重构：EMSeq 相关页职责重排，统一主线状态、`host_mito.CG.cov` 来源与配置必填项说明
- 新增 `saturation` stage：`... -> call -> saturation -> summary`，输出饱和度曲线与汇总表至 `work/<sample>/qc/saturation/`
- `summary/sample_summary.tsv` 新增 `saturation_rate` 字段：从 `qc/saturation/saturation_summary.tsv` 读取；缺失时输出 `NA`
- `workflow/*.json` 新增 `saturation_script` / `saturation_reads_threshold` 配置项及 `slurm.saturation` 资源配置

## 1.3.0

- `call`：按 batch 流式追加写出 `.CG.cov`，控制完整样本内存占用
- `call`：R1/R2 trimming 方向语义与 `mbias` 严格一致，修剪后 CpG coverage 与甲基化统计对齐
- `summary`：重构 sample-level reads 指标，区分 raw / barcoded / host_mapped / host_valid 层级，补齐 spike-in mapped reads
- `summary`：收紧 host mapped reads 与 valid reads 计数口径，排除低置信度或多重比对 reads
- `summary`：per-spot CpG 统计与 heatmap 基于新 per-spot coverage 分布
- 新增 Slurm workflow profile，支持集群环境端到端最小闭环
- 重构首次用户文档阅读路径

## 1.2.0

- `all + slurm + submit` 改为客户端一次性提交完整依赖 DAG，不再依赖计算节点 nested `sbatch`
- `demux_extract_bc` / `align` 在 `slurm` 生成时基于 `number_of_split_parts` 预推导 chunk 作业
- `mbias` / `call` spike 侧脚本以 `spike_in_index` 为权威来源；host `call` 保留运行时输入发现
- 更新 `README.md` 与 `TEST.md`，补充无 nested submit 提交方式与回归检查

## 1.1.3

- 文档明确 `call` 用户侧主交付物，强调 `coverage/*.CG.cov` 与 summary 结果同为一等输出
- `split` stage Slurm 模式固定 `05_split_bams.sbatch -> afterok -> 05_split_sort.sbatch`
- `workflow/dbit_taps_test.json` 调整 Slurm 默认资源请求，`align` CPU 更新为 8
- `scripts/make_cmd.py` 的 `slurm + stage=all` 改为运行时逐 stage 生成并提交下游作业
- `run.sbatch` 与阶段 launcher 基于上游 `final_slurm_dependency` 用 `afterok` 串联
- 补充 README、TEST 与 `doc/commands.md`，说明新 Slurm launcher 行为

## 1.1.2

- `scripts/make_cmd.py` 默认优先解析 `pixi` 环境内工具路径，减少对外部 `PATH` 与 `module load` 依赖
- `fastp_split` Slurm 脚本不再默认注入 `module load fastp`
- `split` stage Slurm 模式新增 `05_split_submit.sh` 作为提交入口

## 1.1.1

- `scripts/extract_bc.py` 默认进度输出间隔从 `1s` 调整为 `60s`
- `workflow/dbit_taps_test.json` 默认同时生成 host 与 spike-in 的 `mbias` 与 `call` 结果

## 1.1.0

- `summary` stage 新增 3 个 heatmap 输出：`reads_heatmap.png`、`cpg_site_count_heatmap.png`、`mean_methylation_heatmap.png`
- 补充 `summary` 相关文档与测试说明

## 1.0.0

- 固定主流程：`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> summary`
- 支持 `local` 与 `Slurm` 运行方式，通过 `pixi` 管理环境
- 交付 `per_spot_summary.tsv`、`sample_summary.tsv` 与 `mbias` 结果作为 1.0 用户侧最小闭环
