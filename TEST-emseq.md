# TEST-EMSeq (Maintainer)

本页是 `EMSeq` 独立入口的维护者回归文档。当前覆盖 `fastp_split -> demux_extract_bc -> align -> pool -> split -> call`，不包含 `saturation -> summary` 的后续链路。用户入口说明请看 `doc/emseq.md`。

## 当前范围

- 入口脚本：`scripts-emseq/make_cmd.py`
- 当前支持 stage：`fastp_split`、`demux_extract_bc`、`align`、`pool`、`split`、`call`
- 单步执行实现：`scripts/fastp_split.py`、`scripts-emseq/extract_bc.py`、`scripts-emseq/aligner.py`、`scripts/pool.py`、`scripts/split_bams.py`、`scripts/bam_sort_parallel.py`
- 单步执行实现：`scripts-emseq/call.py`
- 样例配置：`workflow/dbit_emseq_test.json`

## 最小验收

提交前建议至少执行以下 3 组检查：

1. CLI 可用。
2. 六个 stage 在 `local` 与 `slurm` 下都能正确 dry-run。
3. 用样例配置完成一次 `fastp_split -> demux_extract_bc -> align -> pool -> split -> call` 的本地真实生成。

CLI：

```bash
pixi run python scripts-emseq/make_cmd.py --help
```

`local` dry-run：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage fastp_split \
  --runner local \
  --dry-run

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage demux_extract_bc \
  --runner local \
  --dry-run

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage align \
  --runner local \
  --dry-run

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage pool \
  --runner local \
  --dry-run

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage split \
  --runner local \
  --dry-run

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage call \
  --runner local \
  --dry-run
```

`slurm` dry-run：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage fastp_split \
  --runner slurm \
  --dry-run

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage demux_extract_bc \
  --runner slurm \
  --dry-run

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage align \
  --runner slurm \
  --dry-run

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage pool \
  --runner slurm \
  --dry-run

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage split \
  --runner slurm \
  --dry-run

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage call \
  --runner slurm \
  --dry-run
```

通过标准：

- 命令正常展开
- `fastp_split` / `demux_extract_bc` / `align` / `pool` / `split` / `call` 分别调用对应脚本（含 `scripts-emseq/call.py`）
- 输出目录分别指向 `shard_fastq`、`demux`、`align_shards`、`pooled`、`split_bams`、`pileup`、`coverage`
- 无参数缺失、无路径解析错误
- `slurm` dry-run 中工具路径应解析到当前 `pixi` 环境，或显式使用用户传入的 `--*-bin`

## 真实生成 Smoke

```bash
rm -rf work/test-DNAme-EMSeq

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage fastp_split \
  --runner local \
  --submit

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage demux_extract_bc \
  --runner local \
  --submit

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage align \
  --runner local \
  --submit

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage pool \
  --runner local \
  --submit

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage split \
  --runner local \
  --submit

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage call \
  --runner local \
  --submit
```

## 输出检查点

Smoke 通过标准：

- `work/<sample>/shard_fastq/*.R1.fq.gz`
- `work/<sample>/shard_fastq/*.R2.fq.gz`
- `work/<sample>/shard_fastq/fastp.html`
- `work/<sample>/shard_fastq/fastp.json`
- `work/<sample>/demux/*.R1.demux.fq.gz`
- `work/<sample>/demux/*.R2.demux.fq.gz`
- `work/<sample>/demux/*.R1.spike-in.fq.gz`
- `work/<sample>/demux/*.R2.spike-in.fq.gz`
- `work/<sample>/demux/*.stats.json`
- `work/<sample>/align_shards/*.cb.bam`
- 若配置了 `spike_in_index`：`work/<sample>/align_shards/*.<spike_name>.bam`
- `work/<sample>/pooled/pooled.byCB.bam`
- 若配置了 `spike_in_index`：`work/<sample>/pooled/pooled.<spike_name>.sorted.bam`
- `work/<sample>/split_bams/per_spot_read_counts.tsv`
- `work/<sample>/split_bams/**/*.sorted.bam`
- `work/<sample>/pileup/**/*.vcf.gz`
- `work/<sample>/pileup/**/*.vcf.gz.tbi`
- `work/<sample>/coverage/host/**/*.CG.cov`
- `work/<sample>/coverage/host_mito.CG.cov`
- 若配置了 `spike_in_index`：`work/<sample>/coverage/<spike_name>.CG.cov`

关键行为检查：

- chunk FASTQ 命名应保持与现有 `fastp_split` 契约一致
- `align` 应使用 `demux/*.R1.demux.fq.gz` / `*.R2.demux.fq.gz` 作为 host 输入
- 若配置了 `spike_in_index`，`align` 应先处理 `*.spike-in.fq.gz`，再处理 host 输入
- 若 `spike_in_index` 为空，`pool` 应仅处理 host（不生成 spike/all 相关作业）
- `coverage/host_mito.CG.cov` 应由 `coverage/host/**/*.CG.cov` 中的线粒体 contig 位点汇总得到，默认仅汇总 `chrM`
- `coverage/host/**/*.CG.cov` 不应再包含线粒体 contig 位点，默认不应包含 `chrM`
- EMSeq 入口当前不应暴露 `all` 或其他非 `fastp_split` / `demux_extract_bc` / `align` / `pool` / `split` / `call` stage

## 数据说明

维护者在回归前，建议先确认 `workflow/dbit_emseq_test.json` 里的 `r1` / `r2` 路径与测试目录中的实际文件一致。  
当前仓库中的 `test-DNAme-EMSeq` 目录可见原始数据文件名为：

- `data/raw/test-DNAme-EMSeq/test_ME11_50um_DNAm_1.fq.gz`
- `data/raw/test-DNAme-EMSeq/test_ME11_50um_DNAm_2.fq.gz`
