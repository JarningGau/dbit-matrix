# EMSeq Guide

EMSeq 独立入口用户说明。Stage 输入输出契约见 `doc/stages.md`；结果解读见 `doc/outputs.md`；维护者回归见 `TEST-emseq.md`。

## 支持范围

- **状态**：主线已闭环，接口与产物契约已固定，性能参数可持续调优
- **入口**：`scripts-emseq/make_cmd.py`
- **主流程**：`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> saturation -> summary`
- **`--stage all`**：生成 `work/<sample>/commands/run.sh`（`local`）或 `run.sbatch`（`slurm`），按序串联；与 `mbias_mode` / `call_mode` 的 `all` 取值含义不同
- **`mbias`**：`scripts-emseq/mbias.py`（bisulfite 风格；参考 asTair TOP/BOT，仅统计参考 `CG`；与 TAPS 规则不同）
- **`summary` / `saturation`**：复用 `scripts/summary.py`、`scripts/saturation.py`，产物路径与 TAPS 一致

## 与 TAPS 的主要差异

- 编排入口：`scripts-emseq/make_cmd.py`（非 `scripts/make_cmd.py`）
- `demux_extract_bc`：`scripts/extract_bc.py`
- demux 使用 `linker_bc` 定位 barcode（barcode 长度来自 whitelist），再从 `bc1_end` 之后搜索 `insert_left` 并 trim
- `align`：`scripts-emseq/aligner.py`（`biscuit align`），输出契约与 TAPS 一致（host：`*.cb.bam`；spike：`*.<spike_name>.bam`）
- `call`：默认将 `call_mito_chromosomes`（默认 `chrM`）从 per-spot `coverage/host/**/*.CG.cov` 剔除；`coverage/host_mito.CG.cov` 优先由 `qc/mbias/host.subsampled.sorted.bam` 经 pileup 生成，否则从 per-spot coverage 汇总线粒体位点

## 最小配置（全流程）

从 `workflow/dbit_emseq_test.json` 复制，按 stage 校验必填项：

| 阶段 | 必填字段 |
|------|----------|
| 通用 | `sample_id`、`r1`、`r2`、`work_root` |
| `fastp_split` | `fastp_threads`、`number_of_split_parts`、`fastp_bin` |
| `demux_extract_bc` | `barcode1_whitelist`、`barcode2_whitelist`、`linker_bc`、`insert_left` |
| `align` | `biscuit_reference`、`biscuit_threads`、`biscuit_batch_size`、`biscuit_bin`、`sinto_bin`、`samtools_bin` |
| `pool` | `samtools_threads`、`host_sort_mem`（按需） |
| `split` | **`split_barcodes`**（及 `split_cb_tag`、`split_smoke` 按需） |
| `mbias` | `mbias_script`（默认 `scripts-emseq/mbias.py`）、`mbias_mode`；`all`/`host` 需 `call_reference_file`；`all`/`spike` 需 `spike_in_index` |
| `call` | **`call_reference_file`**、**`call_jobs`**、`call_mode`、`call_mito_chromosomes` |
| `saturation` / `summary` | `saturation_script`、`saturation_reads_threshold`、`summary_script`；`spike_in_index` 供 summary 传 spike 名 |

### 可调参数

- **demux**：`linker_edit_distance`、`barcode_hamming_distance`、`gzip_level`
- **align**：`spike_in_index`（对象或 `NAME=INDEX` 列表）
- **mbias**：`mbias_host_subsample_fraction`、`mbias_max_cycle`、`mbias_min_mapping_quality`；Slurm：`slurm.mbias.host` / `slurm.mbias.spike`
- **call**：`call_host_threads`、`call_spike_threads`、`call_host_subsample_fraction`、`call_host_subsample_seed`（`host_mito` BAM 回退用）；Slurm：`slurm.call.host` / `slurm.call.spike`
- **saturation**：`saturation_reads_threshold`；Slurm：`slurm.saturation`
- **summary**：Slurm：`slurm.summary`
- **`all` + `slurm`**：`run.sbatch` SBATCH 头默认取自 `slurm.summary`

## 常用命令

全流程（建议先 dry-run）：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage all \
  --runner local \
  --dry-run
```

Slurm dry-run：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage all \
  --runner slurm \
  --dry-run
```

单 stage：`--stage` 改为 `fastp_split`、`demux_extract_bc`、`align`、`pool`、`split`、`mbias`、`call`、`saturation`、`summary`；Slurm 时 `--runner local` 改为 `--runner slurm`。

非 dry-run 且 `--submit` 时按 runner 执行或提交生成的脚本。`split` 在 Slurm 下为两段作业链式依赖，见 `doc/stages.md`。

## 结果检查

- **优先查看**：`summary/sample_summary.tsv`、`summary/per_spot_summary.tsv`、`qc/saturation/saturation_summary.tsv`
- **字段说明与排错**：`doc/outputs.md`
- **路径与契约**：`doc/stages.md`

## 常见问题排查

- `reads` 异常：demux 保留率、`split` 输入、`demux/*.stats.json`
- 无 spike 结果：`spike_in_index` 与 `align`/`pool`/`call_mode`
- `host_mito.CG.cov` 或线粒体相关：是否执行 `mbias`（host）、`call` 是否成功；回退路径见「与 TAPS 的主要差异」
- `saturation_rate` 为 `NA`：是否执行 `saturation`，以及 `coverage/host/**/*.CG.cov` 与 `split_bams/per_spot_read_counts.tsv`
