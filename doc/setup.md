# Setup

## 环境

初始化环境：

```bash
pixi lock
pixi install
```

执行命令统一使用 `pixi run`：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --dry-run
```

如需增加 PyPI 依赖：

```bash
pixi add --pypi <package>
pixi lock
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

pixi run make-cmd-dry-run
```

## 测试数据

- R1: `data/raw/test-DNAme-TAPS/250812_intestine_dbit_taps_R1.fq.gz`
- R2: `data/raw/test-DNAme-TAPS/250812_intestine_dbit_taps_R2.fq.gz`

样例 workflow：

- `workflow/dbit_taps_test.json`

## 文库结构

- R1: `barcodeB-linker2-barcodeA-linker1-Tn5-insert`
- R2: `insert`

## 关键参数

- barcodes: `configs/barcodes_50a.tsv`
- linker1: `GTGGCCGATGTTTCG`
- linker2: `ATCCACGTGCTTGAGAGGCCAGAGCATTCG`
- tn5: `CATCGGCGTACGACTAGATGTGTATAAGAGACAG`
- `linker_edit_distance`: linker 模糊匹配编辑距离上限
- `barcode_hamming_distance`: whitelist 纠错汉明距离上限
- `gzip_level`: demux FASTQ 压缩等级
- `bwa_index`: host 参考索引
- `bwa_threads`: align 线程数
- `bwa_bin / sinto_bin / samtools_bin`: 相关工具路径
- `spike_in_index`: spike-in 参考，支持对象或 `NAME=INDEX` 列表
