# DBiT-Matrix

`DBiT-Matrix` 是面向 `DBiT-DNAme-TAPS` 的 1.0 工作流：从原始 `FASTQ` 出发，完成条码提取、比对、spot 拆分、M-bias 质控、甲基化 calling，并输出 spot 级与 sample 级汇总结果。

1.0 版本只聚焦一条可复现、可运行、可交付结果的主流程，不再在 README 中展开开发里程碑或实验性分支。

## 1.0 支持范围

- 当前协议：`DBiT-DNAme-TAPS`
- 固定主流程：`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> summary`
- 运行方式：支持 `local` 和 `Slurm`
- 环境管理：统一使用 `pixi`
- 配置入口：优先通过 `workflow/*.json` 管理参数
- 安全检查：所有 stage 支持 `--dry-run`

## 你会得到什么

对大多数用户，1.0 版本最重要的交付物是：

- `work/<sample>/summary/per_spot_summary.tsv`
- `work/<sample>/summary/sample_summary.tsv`
- `work/<sample>/qc/mbias/*.mbias.tsv`
- `work/<sample>/qc/mbias/*.mbias.png`

其中：

- `per_spot_summary.tsv` 给出每个 spot 的平均甲基化、CpG 位点数和 reads
- `sample_summary.tsv` 给出样本级 host、mito 和 spike-in 的汇总结果
- `mbias` 结果用于检查末端偏倚，不会自动修改 calling 参数

## 输入要求

运行 1.0 工作流前，通常需要准备：

- 原始双端测序数据：`R1 FASTQ`、`R2 FASTQ`
- barcode 白名单：`barcode1_whitelist`、`barcode2_whitelist`
- host 参考序列：`bwa_index` 与 `call_reference_file`
- 可选 spike-in 参考序列：`spike_in_index`
- 一个 workflow 配置文件，例如 `workflow/dbit_taps_test.json`

建议直接从现有样例配置复制一份，再按你的数据路径和参考基因组修改。

## 快速开始

安装环境：

```bash
pixi lock
pixi install
```

先做一次 dry-run，确认配置和命令展开正确：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner local \
  --dry-run
```

生成本地执行脚本：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner local
```

直接提交本地串行入口：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner local \
  --submit
```

提交 Slurm 工作流：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --runner slurm \
  --submit
```

说明：

- 未显式传 `--stage` 时，优先读取 `workflow/*.json` 中的 `stage`
- 推荐在完整流程中显式使用 `--stage all`
- `all` 会依次展开：`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> summary`

## 主流程说明

| Stage | 作用 | 关键输出 |
| --- | --- | --- |
| `fastp_split` | 质控并按 chunk 切分原始 FASTQ | `shard_fastq/*.fq.gz` |
| `demux_extract_bc` | 提取 barcode，分流 matched 与 spike-in reads，并输出统计 | `demux/*.demux.fq.gz`、`demux/*.spike-in.fq.gz`、`demux/*.stats.json` |
| `align` | 对 host 与 spike-in 参考序列比对 | `align_shards/*.cb.bam`、`align_shards/*.<spike_name>.bam` |
| `pool` | 合并 host 与 spike-in 的分片 BAM | `pooled/pooled.byCB.bam`、`pooled/pooled.<spike_name>.sorted.bam` |
| `split` | 按 `CB:Z:<x>+<y>` 将 host BAM 拆成 spot BAM 并排序 | `split_bams/**/*.sorted.bam` |
| `mbias` | 生成 host/spike-in 的 M-bias QC 结果 | `qc/mbias/*` |
| `call` | 进行 host per-spot、host mito 和 spike-in 甲基化 calling | `coverage/*.CG.cov` |
| `summary` | 汇总为 spot 级与 sample 级结果 | `summary/per_spot_summary.tsv`、`summary/sample_summary.tsv` |

## 配置建议

1. 复制 `workflow/dbit_taps_test.json` 作为你的项目配置。
2. 先只改输入路径、样本名、参考基因组和 barcode 白名单。
3. 先执行 `--dry-run`，再执行真实任务。
4. 首次跑通时，优先使用 `local` 或较小测试数据确认输出契约。

如果你只需要某一个 stage，也可以把 `stage` 写进配置文件，或在命令行显式传 `--stage <name>`。

## 输出目录约定

默认工作目录在 `work/<sample>/` 下，常见结构如下：

```text
work/<sample>/
├── shard_fastq/
├── demux/
├── align_shards/
├── pooled/
├── split_bams/
├── qc/mbias/
├── coverage/
└── summary/
```

这套目录约定也是下游 stage 的输入契约，建议不要手动改名或移动中间产物。

## 常见使用场景

只检查命令是否会正确生成：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage all \
  --dry-run
```

只跑某一个 stage：

```bash
pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage mbias \
  --runner local \
  --dry-run
```

查看 CLI 帮助：

```bash
pixi run python scripts/make_cmd.py --help
```

## 文档导航

- `doc/setup.md`：环境、测试数据、样例配置
- `doc/stages.md`：各 stage 的输入输出契约
- `doc/commands.md`：`make_cmd.py` 的常用命令
- `TEST.md`：最小回归与 smoke 检查

## 说明

- 文档示例统一使用 `pixi run ...`
- 如需新增依赖，请同步更新 `pixi.lock`
- 1.0 README 只保留用户使用路径；实现细节与开发约束不再放在首页
