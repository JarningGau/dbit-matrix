# TEST-EMSeq (Maintainer)

EMSeq 独立入口维护者回归。用户说明见 `doc/emseq.md`；输入输出契约见 `doc/stages.md`。

## 范围

- 入口：`scripts-emseq/make_cmd.py`
- Stages：`fastp_split`、`demux_extract_bc`、`align`、`pool`、`split`、`mbias`、`call`、`saturation`、`summary`，实验性 `aggregate`（不纳入 `--stage all`），以及 `all`（生成 `commands/run.sh` / `run.sbatch` 串联全流程）
- 单步脚本：`scripts/fastp_split.py`、`scripts/extract_bc.py`、`scripts-emseq/aligner.py`、`scripts/pool.py`、`scripts/split_bams.py`、`scripts/bam_sort_parallel.py`、`scripts-emseq/mbias.py`、`scripts-emseq/call.py`、`scripts/saturation.py`、`scripts/summary.py`、`scripts/aggregate.py`（`aggregate` stage）
- 样例配置：`workflow/dbit_emseq_test.json`

## 最小验收

提交前至少完成：

1. CLI `--help` 可用
2. 主线九个 stage 与实验性 `aggregate` 在 `local` 与 `slurm` 下均能 dry-run；`--stage all` 在 `local` 与 `slurm` 下 dry-run（不落盘；`all` 不含 `aggregate`）
3. 用样例配置完成一次全流程本地真实生成（可逐 stage，也可 `--stage all` 生成 `run.sh` 后执行）

### CLI

```bash
pixi run python scripts-emseq/make_cmd.py --help
```

### Dry-run：`local`（主线九个 stage + `aggregate` + `all`）

```bash
CFG=workflow/dbit_emseq_test.json
STAGES="fastp_split demux_extract_bc align pool split mbias call saturation summary aggregate"
for stage in $STAGES; do
  pixi run python scripts-emseq/make_cmd.py \
    --workflow-config "$CFG" \
    --stage "$stage" \
    --runner local \
    --dry-run
done

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config "$CFG" \
  --stage all \
  --runner local \
  --dry-run
```

### Dry-run：`slurm`（主线九个 stage + `aggregate` + `all`）

将上述命令中 `--runner local` 改为 `--runner slurm`。

### 通过标准（dry-run）

- 命令正常展开，无缺失参数、无路径解析错误
- 各 stage 调用对应脚本（含 `scripts-emseq/mbias.py`、`scripts-emseq/call.py`、`scripts/saturation.py`、`scripts/summary.py`）
- 输出目录指向 `shard_fastq`、`demux`、`align_shards`、`pooled`、`split_bams`、`qc/mbias`、`pileup`、`coverage`、`qc/saturation`、`summary` 等预期位置
- `slurm` dry-run 中工具路径解析到当前 `pixi` 环境，或显式传入 `--*-bin`

## 真实生成 Smoke（本地）

```bash
rm -rf work/test-DNAme-EMSeq

CFG=workflow/dbit_emseq_test.json
STAGES="fastp_split demux_extract_bc align pool split mbias call saturation summary"
for stage in $STAGES; do
  pixi run python scripts-emseq/make_cmd.py \
    --workflow-config "$CFG" \
    --stage "$stage" \
    --runner local \
    --submit
done
```

### 输出检查点

Smoke 通过时至少存在（路径相对于 `work/<sample>/`）：

- `shard_fastq/*.R1.fq.gz`、`*.R2.fq.gz`、`fastp.html`、`fastp.json`
- `demux/*.R1.demux.fq.gz`、`*.R2.demux.fq.gz`、`*.R1.spike-in.fq.gz`、`*.R2.spike-in.fq.gz`、`*.stats.json`
- `align_shards/*.cb.bam`；配置 `spike_in_index` 时：`align_shards/*.<spike_name>.bam`
- `pooled/pooled.byCB.bam`；spike：`pooled/pooled.<spike_name>.sorted.bam`
- `split_bams/per_spot_read_counts.tsv`、`split_bams/**/*.sorted.bam`
- `qc/mbias/host.subsampled.sorted.bam`（及 `.bai`，`mbias_mode` 含 `host` 时）、`host.mbias.tsv`、`host.mbias.png`
- `pileup/**/*.vcf.gz`（及 `.tbi`）
- `coverage/host/**/*.CG.cov`、`coverage/host_mito.CG.cov`；spike：`coverage/<spike_name>.CG.cov`
- `qc/saturation/saturation_curve.png`、`saturation_summary.tsv`
- `summary/per_spot_summary.tsv`、`sample_summary.tsv`、`*.heatmap.png`

### 关键行为检查

- chunk FASTQ 命名与 `fastp_split` 契约一致
- `align`：host 输入为 `demux/*.demux.fq.gz`；配置 `spike_in_index` 时先 spike-in 后 host
- `spike_in_index` 为空时 `pool` 仅处理 host
- `mbias` 调用 `scripts-emseq/mbias.py`；样例 `mbias_mode: host` 时 local 生成 `commands/06_mbias.sh`（具体以当前编排编号为准）
- `mbias` 为 asTair TOP/BOT 风格，仅参考 `CG` 位点计数
- 若存在 `qc/mbias/host.subsampled.sorted.bam`，`call` 优先用它生成 `coverage/host_mito.CG.cov`；否则由 per-spot coverage 汇总线粒体位点（默认 `chrM`）
- `coverage/host/**/*.CG.cov` 不含线粒体 contig 位点（默认不含 `chrM`）
- `stage=all` 时 `--dry-run` 不落盘；非 dry-run 时生成 `commands/run.sh`（local）或 `run.sbatch`（slurm，client-side sbatch DAG；`split` 仍为两作业链式依赖）
- `--stage all` 的顺序不变，且不包含 `aggregate`

### `aggregate`（实验性）

```bash
pixi run python scripts/aggregate.py --help

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage aggregate \
  --runner local \
  --dry-run

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage aggregate \
  --runner slurm \
  --dry-run
```

单步 smoke（需在 `call` 之后存在 `coverage/host/**/*.CG.cov`）：

```bash
pixi run python scripts/aggregate.py \
  --work-path work/test-DNAme-EMSeq \
  --dry-run
```

检查点：契约与排序规则见 `doc/stages.md` §11；编排产物为 `commands/11_aggregate.sh` 或 `11_aggregate.sbatch`。

## 数据说明

回归前确认 `workflow/dbit_emseq_test.json` 中 `r1` / `r2` 与仓库内数据一致。示例原始文件：

- `data/raw/test-DNAme-EMSeq/test_ME11_50um_DNAm_1.fq.gz`
- `data/raw/test-DNAme-EMSeq/test_ME11_50um_DNAm_2.fq.gz`
