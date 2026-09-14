# Coordinate-sorted pool, splitless mainline, and EMSeq caller

Date: 2026-09-14  
Status: approved  
Scope: TAPS, TAPS v2, and EMSeq `pool`, `call`, `mbias`, `saturation`, `summary`, and all three `make_cmd` drivers

This spec supersedes the follow-ups listed in
[`2026-09-14-taps-call-alignment-design.md`](2026-09-14-taps-call-alignment-design.md)
(skip `split`, port caller to EMSeq). The earlier spec remains the source of
truth for TAPS classification rules, per-spot `.cov` columns, and interval
parallelism. This spec changes the host BAM contract, removes `split` from
the mainline, and adds an EMSeq caller.

## Goal

Sort and index the pooled host BAM by genomic coordinate, drop `split` from
`--stage all`, and have `call` write `coverage/per_spot_read_counts.tsv` from a
full BAM scan. Keep TAPS calling on the shared core. Add
`scripts-emseq/methy_caller.py` with template EMSeq polarity (retained base =
methylated). Do not add `--assay` to `scripts/methy_caller.py`.

## Non-goals

- Switching primary host outputs to template combined `host.CG.cov`
- Changing `.cov` column layouts or pysam 0-based coordinates
- Deleting `scripts/split_bams.py`
- Cabernet / SmC
- Runtime import of `template/DBiT-spatial-Methylation`
- Applying a second C-to-T transform to the barcode whitelist
- Bit-identical equality with biscuit `vcf2bed` / `mergecg` coverage
- Reworking `aggregate` or `methscan` beyond the new host BAM / counts paths
  they already read through `coverage/host/`

## Current state

TAPS / TAPS v2 host `call` already pileups `pooled/pooled.byCB.bam` in genomic
intervals and still requires `split_bams/per_spot_read_counts.tsv`. Host `pool`
sorts by the `CB` tag and does not index. Spike `pool` already coordinate-sorts
and indexes.

EMSeq `call` still runs biscuit pileup on `split_bams/**/*.sorted.bam` and
writes `pileup/**/*.vcf.gz`. Trim is biscuit `-5` / `-3`
(`call_left_trimming` / `call_right_trimming`).

All three drivers include `split` in `STAGE_SEQUENCE` and `STAGE_CHOICES`.
`saturation` and `summary` read `split_bams/per_spot_read_counts.tsv`.
`mbias` and `summary` read `pooled/pooled.byCB.bam`.

Reference caller:
`template/DBiT-spatial-Methylation/script/steps/python/06.methy_caller.py`.

Biscuit baseline (do not overwrite): `work/test-DNAme-EMSeq`.
FASTQ: `data/raw/test-DNAme-EMSeq`.

## Architecture

New mainline for TAPS, TAPS v2, and EMSeq:

`fastp_split -> demux_extract_bc -> align -> pool -> mbias -> call -> saturation -> summary`

`split` is not in `STAGE_SEQUENCE` or `STAGE_CHOICES`. `scripts/split_bams.py`
stays in the repository and is not wired to any workflow driver.

```text
align_shards/*.cb.bam
        │  pool: cat → coordinate sort → index
        ▼
pooled/pooled.host.sorted.bam[.bai|.csi]
        │
        ├─ mbias
        │
        └─ call host
             1. full BAM scan by whitelist CB
                → coverage/per_spot_read_counts.tsv (all assigned spots)
             2. --smoke: write .cov only for top 16 spots
             3. interval pileup
                → coverage/host/<XX>/<XX>_<YY>.CG.cov
                → optional .CH.cov

pooled/pooled.<spike>.sorted.bam     → call aggregate → coverage/<spike>.*.cov
qc/mbias/host.subsampled.sorted.bam  → call aggregate → coverage/host_mito.*.cov
```

### Components

| Unit | Role |
|------|------|
| `scripts/pool.py` | Host: concatenate, coordinate-sort (no `-t CB`), index. Output `pooled/pooled.host.sorted.bam`. Spike unchanged. |
| `scripts/methy_call_core.py` | Shared: interval pileup, CB → `XX_YY`, full-BAM counts, smoke selection, host shuffle, index checks. Internal polarity flag only; no `--assay` CLI. |
| `scripts/methy_caller.py` | TAPS entry: converted base = methylated. `--mode host\|aggregate`. |
| `scripts-emseq/methy_caller.py` | EMSeq entry: retained base = methylated. Same modes. Imports core; does not import TAPS classify functions as the polarity source. |
| `scripts/call.py` | TAPS / TAPS v2 orchestrator. Host BAM is `pooled.host.sorted.bam`. Does not read `split_bams/`. |
| `scripts-emseq/call.py` | EMSeq orchestrator. Same host/mito/spike layout. No biscuit pileup or VCF. |
| Three `make_cmd.py` | Drop `split`. Host checks: pooled host BAM + index + barcodes. `split_smoke` → `call_smoke`. |
| `saturation` / `summary` / `mbias` | New counts path and new host BAM name. |

CB rules stay those of `scripts/split_bams.py`, not the template spot-ID scheme:

- tag `CB` (configurable, default `CB`)
- `CB:Z:<x>+<y>`, or concatenated `xy` of length `2 * barcode_len`
- whitelist: first column, uppercase, 0-based order
- stem `{x_index:02d}_{y_index:02d}` under `coverage/host/{x_index:02d}/`
- missing tag, unparsable CB, or barcode not in the whitelist: skip, do not raise

Do not apply template whitelist C-to-T. Demux already writes the barcode
sequence that the whitelist uses (EMSeq test configs already point at
`configs/barcodes_50a_CtoT.tsv`).

### Parallelism

Unchanged from the TAPS pooled-caller spec: `call_jobs` is interval workers
inside one chromosome (default 8). `batch_size` default 10_000_000 bp.
Slurm: one host sbatch and one sbatch per spike. Set
`slurm.call.host.cpus_per_task` to at least `call_jobs`.

## Calling semantics

### Shared filters

Same as the approved TAPS spec and the template:

- Keep properly paired primary flags `{99, 147, 83, 163}` only
- Skip deletion, refskip, missing `query_position`, missing sequence
- Pileup: `ignore_overlaps=True`, `ignore_orphans=True`, `stepper="samtools"`,
  `redo_baq=False`, `compute_baq=True`, `truncate=True`
- Missing `YD`: keep. `YD == "u"`: skip
- Use `query_sequence`, not `get_forward_sequence()`
- Defaults: `min_base_quality=30`, `min_mapping_quality=10`, `max_depth=50000`
- R1/R2 left/right trim uses the existing cycle definition in
  `scripts/methy_caller.py`. CG trims both dinucleotide query positions. CH
  trims the cytosine query position; plus-strand CH also trims the neighbor
  query position

EMSeq drops biscuit `-5` / `-3`. It uses the same four trim knobs as TAPS
(`call_trimming.host` / `call_trimming.spike`, each with r1/r2 left/right).
Remove `call_left_trimming` and `call_right_trimming`.

### Polarity

Core takes an internal polarity, set only by each entrypoint:

| Observation | TAPS (`converted_means_methylated=True`) | EMSeq (`False`) |
|-------------|------------------------------------------|-----------------|
| CpG `TG` / `CA` | methylated | unmethylated |
| CpG `CG` | unmethylated | methylated |
| CH plus `T` / minus `A` | methylated | unmethylated |
| CH plus `C` / minus `G` | unmethylated | methylated |

CH site discovery, strand flags, and neighbor rules match the TAPS spec
(plus `{99, 147}`, minus `{83, 163}`; CC allows converted neighbor).
`call_context_mode` remains `cg` | `ch` | `both`, default `cg`. EMSeq gains
CH for the first time; default stays `cg`.

One interval pileup classifies both contexts when `both` is set.

### Counts, roster, and smoke

Host `call` scans the pooled host BAM with `fetch(until_eof=True)` and counts
every record whose CB maps to a whitelist spot. Semantics match
`scripts/split_bams.py` assigned reads, not pileup-overlapping reads.

If `coverage/per_spot_read_counts.tsv` already exists and `--force-recount` is
not set, skip the scan and reuse that file. Otherwise write it with header:

```text
X_index	Y_index	spot	reads
```

Rows are all assigned spots, sorted by `X_index`, `Y_index`. Smoke and
`--spot-roster` never shrink this file.

The `.cov` write-set is exactly one of these, in this order:

1. `--spot-roster PATH` if set. The file is either the counts TSV (`spot`
   column) or one `XX_YY` id per line. `--smoke` is ignored when a roster is
   set.
2. Else if `--smoke` / `call_smoke`: the 16 spots with highest `reads`, ties
   broken by `spot` ascending. Deterministic. Do not reuse
   `split --smoke`'s `random.sample`.
3. Else every spot in the counts TSV.

A write-set spot with no passing methylation observations gets an empty
`.cov`. Do not create files for whitelist spots that never appear in the BAM.

`--spot-roster` is for biscuit comparison against the 16 existing baseline
spots. Default: no roster file.

### Aggregate

Mito and spike use `--mode aggregate` (no barcodes, no counts TSV).

Mito input remains `qc/mbias/host.subsampled.sorted.bam` when present, else
subsample from `pooled/pooled.host.sorted.bam`.

EMSeq no longer writes `pileup/**/*.vcf.gz`. EMSeq no longer builds
`host_mito.CG.cov` by stripping mito rows from per-spot files when the
subsample BAM exists; it uses the same aggregate path as TAPS.

## Outputs and compatibility

Host and aggregate `.cov` paths and columns stay as in
[`docs/developers/contracts.md`](../../developers/contracts.md) `call`.
Positions stay pysam 0-based. `pct` stays
`round(100 * meth / (meth + unmeth), 2)` with two decimal places for TAPS
and EMSeq (biscuit currently emits integer percent; that is not preserved).

Write through a process-specific temp path and `os.replace`. Reruns skip per
final host `.cov` that already exists for a write-set spot. Existing mito or
spike `.cov` is not recomputed.

## Orchestration and config

### `pool`

Host commands:

```text
samtools cat -o pooled/pooled.host.unsorted.bam <shards>
samtools sort -m <host_sort_mem> -@ <threads> -o pooled/pooled.host.sorted.bam pooled/pooled.host.unsorted.bam
samtools index pooled/pooled.host.sorted.bam
rm pooled/pooled.host.unsorted.bam
```

No `-t CB`. Spike stays `pooled/pooled.<spike>.sorted.bam` plus index.

### `call` wrappers

Host requires `--barcodes` (from workflow `split_barcodes`, default
`barcode1_whitelist`). The key name `split_barcodes` stays so existing JSON
does not rename; it is the CB whitelist for `call`, not a split-stage input.

Host requires `pooled/pooled.host.sorted.bam` and a `.bai` or `.csi`.
Host does not require `split_bams/`.

EMSeq `scripts-emseq/call.py` mirrors TAPS `scripts/call.py` (host / mito /
spike, trim nest, context mode, jobs). Delete biscuit pileup, `vcf2bed`,
`mergecg`, tabix, `call_host_threads`, `call_spike_threads`,
`call_left_trimming`, and `call_right_trimming`.

### Drivers

`scripts/make_cmd.py`, `scripts-v2/make_cmd.py`, `scripts-emseq/make_cmd.py`:

- Remove `split` from `STAGE_SEQUENCE` and `STAGE_CHOICES`
- `--stage split` is an error
- Drop split command builders, split Slurm helpers, and `slurm.split` from
  example JSON
- Pass `--barcodes` from `split_barcodes` into `call`
- Add `call_smoke` (bool, default false). e2e configs that set
  `split_smoke: true` become `call_smoke: true`
- Host workdir checks: `pooled.host.sorted.bam`, its index, barcode TSV
  (still skip workdir intermediates for `--stage all`)
- EMSeq: add `call_context_mode` and TAPS-style `call_trimming`; drop the
  biscuit call keys listed above
- Keep `call_jobs`, `call_batch_size`, `call_chromosomes`,
  `call_mito_chromosomes`, `call_sample_size`, `call_min_base_quality`,
  `call_min_mapping_quality`, `call_max_depth`

### Downstream

| Consumer | New input |
|----------|-----------|
| `scripts/saturation.py` | `coverage/per_spot_read_counts.tsv` |
| `scripts/summary.py` and `scripts-v2/summary.py` | same TSV; host BAM `pooled/pooled.host.sorted.bam` |
| `scripts/mbias.py` and `scripts-emseq/mbias.py` | `pooled/pooled.host.sorted.bam` |
| `scripts/host_subsample_bam.py` callers | same host BAM name |

`aggregate` and `methscan` still read `coverage/host/**/*.CG.cov` only.

## Errors and dry-run

Non-zero exit when:

- BAM, BAM index, reference FASTA, or `.fai` is missing
- requested chromosomes are absent from the reference
- `--jobs`, `--batch-size`, or `--max-depth` is not `> 0`
- quality or trim values are negative
- host mode is missing `--barcodes`
- barcode whitelist is empty or has unequal lengths

`--dry-run` prints resolved parameters (mode, BAM, reference, chromosomes,
context, jobs, trims, barcodes path when host, smoke) and writes nothing.

## Documentation (same change)

Update these pages; do not copy the same default twice:

- `AGENTS.md` (mainline without `split`; pool output name; Slurm `split`
  layout removed from required stage layout)
- `docs/developers/contracts.md` (pool host name + index; `split` no longer
  mainline; `call` inputs; saturation/summary counts path)
- `docs/developers/config-reference.md` (`call_smoke`, dropped biscuit call
  keys, dropped `split` fields from required mainline)
- `docs/developers/architecture.md` and `docs/developers/runner-reference.md`
  if they still list `split` in `--stage all`
- `docs/users/taps.md` and `docs/users/emseq.md`
- `docs/maintenance/changelog.md` under `## Current`, with `BREAKING`
- `docs/maintenance/taps-regression.md` and
  `docs/maintenance/emseq-regression.md` (stage order; biscuit compare)

`doc/` is not updated unless a task targets it.

The previous TAPS call spec file stays as historical context. Do not rewrite
it in place; this file is the active contract for the new pool / split /
EMSeq work.

## Verification

Required before calling the change done:

- `pixi run python scripts/make_cmd.py --version` and the v2 / EMSeq
  equivalents
- `--help` on `scripts/pool.py`, `scripts/call.py`, `scripts-emseq/call.py`,
  both `methy_caller.py` files
- TAPS, TAPS v2, and EMSeq `make_cmd --stage all --dry-run` (local and, for
  EMSeq, the existing slurm dry-run). Commands must not mention `split`.
  Host pool command must coordinate-sort and index `pooled.host.sorted.bam`

Unit tests (stdlib `unittest`; do not add pytest or change `pixi.lock`):

- TAPS CG: `TG`/`CA` methylated, `CG` unmethylated, other dinucleotides ignored
- EMSeq CG: inverted
- CH plus/minus flag sets and polarity for both assays
- CH CC neighbor allows C/T (plus) and G/A (minus)
- CH CA/CT neighbor must match exactly
- `YD:u` skipped; missing `YD` kept
- CB `x+y` and concatenated forms match `split_bams.py` indices and `XX_YY`
- counts scan: assigned vs skipped
- smoke: full TSV, `.cov` only for top 16 by `reads` desc, `spot` asc
- `--dry-run` does not write outputs

### Biscuit regression

Do not modify `work/test-DNAme-EMSeq`. That directory is the biscuit
baseline (16 host `.CG.cov` files, `pileup/**/*.vcf.gz`,
`pooled/pooled.byCB.bam` without a host index).

FASTQ source: `data/raw/test-DNAme-EMSeq`. Config starting point:
`workflow/dbit_emseq_test.json`.

Procedure:

1. Create a sibling work directory, for example
   `work/test-DNAme-EMSeq-pooled-call/`.
2. Coordinate-sort and index the existing
   `work/test-DNAme-EMSeq/pooled/pooled.byCB.bam` into
   `pooled/pooled.host.sorted.bam` in the sibling (or copy `align_shards` and
   run new `pool`). Copy or point at the same barcode TSV.
3. Run the new EMSeq host `call` with `--spot-roster` listing the 16 spots
   already present under `work/test-DNAme-EMSeq/coverage/host/` so the
   write-set matches biscuit. Do not use smoke selection for this compare;
   new smoke is top-16 by reads and may differ from the old random 16.
4. Compare with a small script checked into the repo (path:
   `scripts-emseq/compare_cov_to_biscuit.py`): join on
   `spot + chrom + coordinate`, mapping biscuit `pos` to 0-based as
   `pos - 1`. Biscuit lines look like
   `chr10 56238827 56238827 100 1 0`.
5. Do not require row-level equality. Print: shared-site fraction, Pearson
   correlation of `meth` and of `unmeth` on shared sites, mean absolute
   `pct` delta, count of sites only in biscuit, count only in the new caller.
6. Pass bar: among shared sites with `meth + unmeth > 0` on both sides, the
   fraction whose majority call disagrees (`meth > unmeth` on one side and
   not the other) must be `< 0.5`. Record shared-site fraction, correlations,
   and mean `|Δpct|` in `docs/maintenance/changelog.md` under `## Current`.
   Biscuit is not a bit-identical gold standard.

A 500 Gbp benchmark is not required.

## Breaking changes

Record explicitly in the changelog body:

- Host pooled BAM is `pooled/pooled.host.sorted.bam`, coordinate-sorted and
  indexed. `pooled.byCB.bam` is no longer written.
- `split` is removed from `--stage all` and from `STAGE_CHOICES`.
- `saturation` / `summary` read `coverage/per_spot_read_counts.tsv`, written
  by `call`.
- `split_smoke` is replaced by `call_smoke` (deterministic top 16).
- EMSeq host/mito/spike calling no longer uses biscuit pileup or VCF
  intermediates. EMSeq CG percentages and site sets will move.
- EMSeq trim keys `call_left_trimming` / `call_right_trimming` are removed;
  use `call_trimming` like TAPS.
- EMSeq `.CG.cov` positions are pysam 0-based, not biscuit's transformed
  coordinate. `pct` has two decimal places.

## Out of scope follow-ups

- Emit template-style combined `host.CG.cov` as the primary artifact
- Re-wire `scripts/split_bams.py` as an optional `make_cmd` stage
- Cabernet / SmC polarity in this repository
- Keep biscuit pileup as a hidden compare mode inside `scripts-emseq/call.py`
  (compare uses the frozen `work/test-DNAme-EMSeq` tree instead)
