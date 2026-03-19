# TEST-EMSeq (Maintainer)

本页是 `EMSeq` 独立入口当前 MVP 的维护者回归文档。  
当前覆盖 `fastp_split -> demux_extract_bc` 两个 stage，不包含 `align -> pool -> split -> mbias -> call` 的后续链路。

如果你是首次使用者，请先阅读：`README.md -> doc/config.md -> doc/stages.md`。

## 当前范围

- 入口脚本：`scripts-emseq/make_cmd.py`
- 当前支持 stage：`fastp_split`、`demux_extract_bc`
- 单步执行实现：`scripts/fastp_split.py`、`scripts-emseq/extract_bc.py`
- 样例配置：`workflow/dbit_emseq_test.json`

## 最小验收

建议在提交前至少执行以下几步：

1. 检查 CLI 是否可用。
2. 检查 `fastp_split` 在 `local` 下能否正确 dry-run。
3. 检查 `fastp_split` 在 `slurm` 下能否正确 dry-run。
4. 检查 `demux_extract_bc` 在 `local` 下能否正确 dry-run。
5. 检查 `demux_extract_bc` 在 `slurm` 下能否正确 dry-run。

CLI 检查：

```bash
pixi run python scripts-emseq/make_cmd.py --help
```

`fastp_split` `local` dry-run：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage fastp_split \
  --runner local \
  --dry-run
```

`fastp_split` `slurm` dry-run：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage fastp_split \
  --runner slurm \
  --dry-run
```

通过标准：

- 命令正常展开
- stage 固定为 `fastp_split` 或显式传入 `demux_extract_bc`
- `fastp_split` 生成命令显式调用 `scripts/fastp_split.py`
- `demux_extract_bc` 生成命令显式调用 `scripts-emseq/extract_bc.py`
- `fastp_split` 输出目录指向 `work/<sample>/shard_fastq`
- `demux_extract_bc` 输出目录指向 `work/<sample>/demux`
- 无参数缺失、无路径解析错误
- `slurm` dry-run 中 `fastp` 默认应展开为当前 `pixi` 环境下的可执行路径，或显式使用用户传入的 `--fastp-bin`

`demux_extract_bc` `local` dry-run：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage demux_extract_bc \
  --runner local \
  --dry-run
```

`demux_extract_bc` `slurm` dry-run：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage demux_extract_bc \
  --runner slurm \
  --dry-run
```

## 真实生成回归

```bash
rm -rf work/test-DNAme-EMSeq

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage fastp_split \
  --runner local
```

通过标准：

- 会写出 `work/test-DNAme-EMSeq/commands/01_fastp_split.sh`
- 脚本内容只负责执行 `fastp_split`
- 不要求预先存在 `work/test-DNAme-EMSeq/shard_fastq`

`slurm` 真实生成回归：

```bash
rm -rf work/test-DNAme-EMSeq

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage fastp_split \
  --runner slurm
```

通过标准：

- 会写出 `work/test-DNAme-EMSeq/commands/01_fastp_split.sbatch`
- `sbatch` 文件只包含 `fastp_split` 所需资源配置
- 日志路径位于 `work/test-DNAme-EMSeq/logs/`

`demux_extract_bc` 真实生成回归（可选）：

```bash
rm -rf work/test-DNAme-EMSeq

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage fastp_split \
  --runner local

pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage demux_extract_bc \
  --runner local
```

通过标准：

- `fastp_split` 执行结束后，`work/test-DNAme-EMSeq/shard_fastq/` 下存在 `*.R1.fq.gz`、`*.R2.fq.gz`
- `demux_extract_bc` 执行结束后，`work/test-DNAme-EMSeq/demux/` 下存在 `*.R1.demux.fq.gz`、`*.R2.demux.fq.gz`、`*.R1.spike-in.fq.gz`、`*.R2.spike-in.fq.gz`、`*.stats.json`

## 可选：直接检查单步脚本

如果需要把问题定位到执行层，而不是编排层，可以直接调用 `fastp_split` 单步脚本：

```bash
pixi run python scripts/fastp_split.py \
  --r1 data/raw/test-DNAme-EMSeq/test_ME11_50um_DNAm_1.fq.gz \
  --r2 data/raw/test-DNAme-EMSeq/test_ME11_50um_DNAm_2.fq.gz \
  --work-path work/test-DNAme-EMSeq.manual \
  --fastp-threads 2 \
  --number-of-split-parts 2 \
  --dry-run
```

如需直接检查 EMSeq 的 barcode 提取脚本，可执行：

```bash
pixi run python scripts-emseq/extract_bc.py \
  work/test-DNAme-EMSeq/shard_fastq/0001.R1.fq.gz \
  work/test-DNAme-EMSeq/shard_fastq/0001.R2.fq.gz \
  -b1 configs/barcodes_50a_CtoT.tsv \
  -b2 configs/barcodes_50a_CtoT.tsv \
  -o work/test-DNAme-EMSeq.manual/demux/0001 \
  --linker1 TAAGTGTTGGTTTTTTGTATTT \
  --linker2 ATTTATGTGTTTGAGAGGTTAGAGTATTTG \
  --tn5 TATTGGTGTATGATTAGATGTGTATAAGAGATAG \
  --linker-edit-distance 1 \
  --barcode-hamming-distance 0
```

通过标准：

- 输出目录为 `work/test-DNAme-EMSeq.manual/shard_fastq`
- 命令中包含 `--split 2`
- 会生成 `fastp.html` 与 `fastp.json` 的输出路径

## 输出检查点

真实运行 `fastp_split` 后，重点检查：

- `work/<sample>/shard_fastq/*.R1.fq.gz`
- `work/<sample>/shard_fastq/*.R2.fq.gz`
- `work/<sample>/shard_fastq/fastp.html`
- `work/<sample>/shard_fastq/fastp.json`

如同时运行 `demux_extract_bc`，补充检查：

- `work/<sample>/demux/*.R1.demux.fq.gz`
- `work/<sample>/demux/*.R2.demux.fq.gz`
- `work/<sample>/demux/*.R1.spike-in.fq.gz`
- `work/<sample>/demux/*.R2.spike-in.fq.gz`
- `work/<sample>/demux/*.stats.json`

关键行为检查：

- chunk FASTQ 命名应保持与现有 `fastp_split` 契约一致
- `fastp.html` / `fastp.json` 应位于 `shard_fastq/`
- EMSeq 入口当前不应暴露 `all` 或其他非 `fastp_split` / `demux_extract_bc` stage

## 数据说明

维护者在回归前，建议先确认 `workflow/dbit_emseq_test.json` 里的 `r1` / `r2` 路径与测试目录中的实际文件一致。  
当前仓库中的 `test-DNAme-EMSeq` 目录可见原始数据文件名为：

- `data/raw/test-DNAme-EMSeq/test_ME11_50um_DNAm_1.fq.gz`
- `data/raw/test-DNAme-EMSeq/test_ME11_50um_DNAm_2.fq.gz`
