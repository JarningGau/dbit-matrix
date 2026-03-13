# Version Log

本页记录面向用户可见的版本变化。`README.md` 只保留当前版本号，具体变化统一写在这里。

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
