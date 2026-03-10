pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage demux_extract_bc \
  --runner local \
  --dry-run

pixi run python scripts/extract_bc.py \
  work/test-DNAme-TAPS/shard_fastq/0001.R1.fq.gz \
  work/test-DNAme-TAPS/shard_fastq/0001.R2.fq.gz \
  -b1 configs/barcodes_50a.tsv \
  -b2 configs/barcodes_50a.tsv \
  -o work/test-DNAme-TAPS/demux/0001.test \
  --linker-edit-distance 1 \
  --barcode-hamming-distance 1 \
  --gzip-level 1

[extract_bc] kept=346137/487575 spike_in=141438 avg_speed=41021.6 reads/s

# 可选：严格基线（仅精确匹配）用于速度/结果对照
pixi run python scripts/extract_bc.py \
  work/test-DNAme-TAPS/shard_fastq/0001.R1.fq.gz \
  work/test-DNAme-TAPS/shard_fastq/0001.R2.fq.gz \
  -b1 configs/barcodes_50a.tsv \
  -b2 configs/barcodes_50a.tsv \
  -o work/test-DNAme-TAPS/demux/0001.strict \
  --linker-edit-distance 0 \
  --barcode-hamming-distance 0 \
  --gzip-level 1

[extract_bc] kept=275554/487575 spike_in=212021 avg_speed=64884.6 reads/s

# align CLI 帮助与命令生成回归
pixi run python scripts/align.py --help

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage align \
  --runner local \
  --dry-run

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage align \
  --runner slurm \
  --dry-run

# 可选：验证多 spike-in 顺序（先 spike-in 再 host）
pixi run python scripts/align.py \
  --work-path work/test-DNAme-TAPS \
  --chunk 0001 \
  --bwa-index /mnt/e/LiLab_HL/resource/bwa/mm10/genome.fa \
  --spike-in-index lambda=/mnt/e/LiLab_HL/resource/bwa/lambda/genome.fa \
  --spike-in-index puc19=/mnt/e/LiLab_HL/resource/bwa/puc19/genome.fa \
  --dry-run

# pool CLI 帮助与命令生成回归
pixi run python scripts/pool.py --help

pixi run python scripts/pool.py \
  --work-path work/test-DNAme-TAPS \
  --samtools-bin samtools \
  --spike-in-name lambda \
  --spike-in-name puc19 \
  --dry-run

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage pool \
  --runner local \
  --dry-run

pixi run python scripts/make_cmd.py \
  --workflow-config workflow/dbit_taps_test.json \
  --stage pool \
  --runner slurm \
  --dry-run