# Version Log

本页记录面向用户可见的版本变化。`README.md` 只保留当前版本号，具体变化统一写在这里。
 
## 1.4.0

- 主流程新增 `saturation` stage：`... -> call -> saturation -> summary`，在 `work/<sample>/qc/saturation/` 下输出饱和度曲线图与汇总表
- `summary/sample_summary.tsv` 新增 `saturation_rate` 字段：从 `qc/saturation/saturation_summary.tsv` 读取；缺失时输出 `NA`
- `workflow/*.json` 新增 `saturation_script` / `saturation_reads_threshold` 配置项，并补齐对应的 `slurm.saturation` 资源配置块

## 1.3.0

- `call` 阶段按 batch 流式追加写出 `.CG.cov`，在完整样本上保持甲基化调用的内存占用在可预测范围内
- `call` 阶段对 R1/R2 trimming 的方向语义与 `mbias` 保持严格一致，修剪后 CpG 的 coverage 与甲基化统计与可视化结果对齐
- `summary` 阶段重构 sample-level reads 指标：统一以 reads 为主度量，区分 raw / barcoded / host_mapped / host_valid 等不同层级，并补齐 spike-in mapped reads
- `summary` 阶段收紧 host mapped reads 与 valid reads 的计数口径，避免将低置信度或多重比对的 reads 计入可用 reads
- `summary` 阶段的 per-spot CpG 统计与 heatmap 基于新的 per-spot coverage 分布，帮助更直观地理解 spot 级 CpG 丰度与甲基化水平
- 新增适配 Slurm 的 workflow profile，使集群环境可以直接从样例配置启动端到端最小闭环
- 重构首次用户文档阅读路径，从 FASTQ 到 summary 保持一条低摩擦主线
 
## 1.2.0

- `all + slurm + submit` 改为由客户端一次性提交完整依赖 DAG，不再依赖计算节点内的 nested `sbatch` 提交
- `demux_extract_bc` / `align` 在 `slurm` 生成时改为基于 `number_of_split_parts` 预推导 chunk 作业，降低对上游中间文件存在性的耦合
- `mbias` / `call` 的 spike 侧脚本改为以 `spike_in_index` 为权威来源；host `call` 脚本改为保留运行时输入发现
- 更新 `README.md` 与 `TEST.md`，补充无 nested submit 的提交方式与回归检查

## 1.1.3

- README、`doc/stages.md` 与 `doc/commands.md` 进一步明确 `call` 的用户侧主交付物，强调 `coverage/*.CG.cov` 与 summary 结果一样属于一等输出
- `split` stage 在单独使用 `Slurm` 时继续固定保持 `05_split_bams.sbatch -> afterok -> 05_split_sort.sbatch`，减少手动提交和直接提交之间的行为差异
- `workflow/dbit_taps_test.json` 调整了样例 Slurm 默认资源请求，使共享集群上的最小闭环验证更容易排队；当前样例配置中的 `align` CPU 请求同步更新为 8
- `scripts/make_cmd.py` 的 `slurm + stage=all` 改为运行时逐 stage 生成并提交下游作业，允许在空 `work/` 目录下直接启动完整流程
- `run.sbatch` 与阶段 launcher 会基于上游 stage 返回的 `final_slurm_dependency` 用 `afterok` 串联，避免在上游输出尚不存在时提前展开下游脚本
- 补充 README、TEST 和 `doc/commands.md`，明确新的 Slurm launcher 行为，并提供提交依赖图便于排查和理解

## 1.1.2

- `scripts/make_cmd.py` 默认优先解析当前 `pixi` 环境内的 `fastp`、`bwa`、`sinto`、`samtools` 可执行路径，减少 `local` 与 `Slurm` 运行时对外部 `PATH` 和 `module load` 的依赖
- `fastp_split` 的 Slurm 脚本不再默认注入 `module load fastp`，并补充对应的 README/TEST 说明
- `split` stage 在单独使用 `Slurm` 时，新增 `05_split_submit.sh` 作为提交入口，固定用 `afterok` 串联 `05_split_bams.sbatch` 与 `05_split_sort.sbatch`
- 补充 `split`/`sort` 阶段文档，明确两段脚本的职责边界与依赖关系

## 1.1.1

- 将 `scripts/extract_bc.py` 的默认进度输出间隔从 `1s` 调整为 `60s`，减少长任务运行时的日志噪声
- 更新 `workflow/dbit_taps_test.json`，默认同时生成 host 与 spike-in 的 `mbias` 和 `call` 结果，便于做端到端最小闭环验收

## 1.1.0

- `summary` stage 新增 3 个基于 `summary/per_spot_summary.tsv` 的 heatmap 输出：
  `reads_heatmap.png`、`cpg_site_count_heatmap.png`、`mean_methylation_heatmap.png`
- 补充 `summary` 相关文档与测试说明，明确 heatmap 的输入来源和输出契约

## 1.0.0

- 固定主流程为 `fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> summary`
- 提供 `local` 和 `Slurm` 两种运行方式，并统一通过 `pixi` 管理环境
- 交付稳定的 `per_spot_summary.tsv`、`sample_summary.tsv` 和 `mbias` 结果，作为 1.0 用户侧最小闭环
