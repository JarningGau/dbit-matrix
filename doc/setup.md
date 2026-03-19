# Setup

本页面向首次用户，聚焦环境准备、测试数据和最小输入材料。配置字段说明请看 `doc/config.md`，运行命令请看 `doc/commands.md`。如果你使用 `EMSeq` 独立入口，请改看 `doc/emseq.md`。

## 环境

初始化环境：

```bash
pixi install
```

执行命令统一使用 `pixi run`：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner local \
  --dry-run
```

## 环境验收

```bash
pixi install

pixi run which python
pixi run python --version
pixi run which fastp
pixi run fastp --version
pixi run which bwa
pixi run bwa 2>&1 | head -n 1
pixi run which sinto
pixi run sinto --help >/dev/null
pixi run which samtools
pixi run samtools --version | head -n 1
```

## 你需要准备什么

- 双端测序数据：`R1 FASTQ`、`R2 FASTQ`
- barcode 白名单：`barcode1_whitelist`、`barcode2_whitelist`
- host 参考：`bwa_index`、`call_reference_file`
- 可选 spike-in 参考：`spike_in_index`
- 一份 workflow 配置文件，建议从 `workflow/dbit_taps_test.json` 复制

## 测试数据

- R1：`data/raw/test-DNAme-TAPS/250812_intestine_dbit_taps_R1.fq.gz`
- R2：`data/raw/test-DNAme-TAPS/250812_intestine_dbit_taps_R2.fq.gz`
- 样例 workflow：`workflow/dbit_taps_test.json`

## 文库结构

TAPS 当前假定：

- R1：`barcodeB-linker2-barcodeA-linker1-Tn5-insert`
- R2：`insert`

EMSeq 的文库结构和最小配置说明见 `doc/emseq.md`。
