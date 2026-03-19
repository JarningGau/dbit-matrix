# TEST-EMSeq (Maintainer)

本页是 `EMSeq` 独立入口当前 MVP 的维护者回归文档。  
当前只覆盖 `fastp_split` stage，不包含 `demux_extract_bc -> align -> pool -> split -> mbias -> call` 的后续链路。

如果你是首次使用者，请先阅读：`README.md -> doc/config.md -> doc/stages.md`。

## 当前范围

- 入口脚本：`scripts-emseq/make_cmd.py`
- 当前唯一 stage：`fastp_split`
- 共享执行实现：`scripts/fastp_split.py`
- 样例配置：`workflow/dbit_emseq_test.json`

## 最小验收

建议在提交前至少执行以下 4 步：

1. 检查 CLI 是否可用。
2. 检查 `fastp_split` 在 `local` 下能否正确 dry-run。
3. 检查 `fastp_split` 在 `slurm` 下能否正确 dry-run。
4. 检查 `fastp_split` 在非 dry-run 模式下能否正确生成命令脚本。

CLI 检查：

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
```

`slurm` dry-run：

```bash
pixi run python scripts-emseq/make_cmd.py \
  --workflow-config workflow/dbit_emseq_test.json \
  --stage fastp_split \
  --runner slurm \
  --dry-run
```

通过标准：

- 命令正常展开
- stage 固定为 `fastp_split`
- 生成命令显式调用 `scripts/fastp_split.py`
- 输出目录指向 `work/<sample>/shard_fastq`
- 无参数缺失、无路径解析错误
- `slurm` dry-run 中 `fastp` 默认应展开为当前 `pixi` 环境下的可执行路径，或显式使用用户传入的 `--fastp-bin`

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
- 不会生成 `demux_extract_bc` 及后续 stage 的命令文件

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

## 可选：直接检查单步脚本

如果需要把问题定位到执行层，而不是编排层，可以直接调用共享的 `fastp_split` 单步脚本：

```bash
pixi run python scripts/fastp_split.py \
  --r1 data/raw/test-DNAme-EMSeq/test_ME11_50um_DNAm_1.fq.gz \
  --r2 data/raw/test-DNAme-EMSeq/test_ME11_50um_DNAm_2.fq.gz \
  --work-path work/test-DNAme-EMSeq.manual \
  --fastp-threads 2 \
  --number-of-split-parts 2 \
  --dry-run
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

关键行为检查：

- chunk FASTQ 命名应保持与现有 `fastp_split` 契约一致
- `fastp.html` / `fastp.json` 应位于 `shard_fastq/`
- EMSeq 入口当前不应暴露 `all` 或其他非 `fastp_split` stage

## 数据说明

维护者在回归前，建议先确认 `workflow/dbit_emseq_test.json` 里的 `r1` / `r2` 路径与测试目录中的实际文件一致。  
当前仓库中的 `test-DNAme-EMSeq` 目录可见原始数据文件名为：

- `data/raw/test-DNAme-EMSeq/test_ME11_50um_DNAm_1.fq.gz`
- `data/raw/test-DNAme-EMSeq/test_ME11_50um_DNAm_2.fq.gz`
