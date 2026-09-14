# TAPS methylation calling alignment with DBiT-spatial-Methylation

Date: 2026-09-14  
Status: draft, pending user review  
Scope: TAPS and TAPS v2 `call` only

## Goal

Replace the TAPS / TAPS v2 Python methylation callers with a template-aligned
implementation: one pileup over `pooled/pooled.byCB.bam`, interval parallelism,
and the template TAPS classification rules. Keep the current per-spot and
aggregate `.cov` paths so `saturation`, `summary`, `aggregate`, and `methscan`
do not change in this iteration.

## Non-goals

- EMSeq / biscuit pileup (`scripts-emseq/call.py`)
- Removing `split` from `--stage all`
- Changing saturation, summary, aggregate, or methscan
- Changing `.cov` column layouts or switching to template combined
  `host.CG.cov` files
- Introducing an `--assay` switch (this caller is TAPS / TAPS v2 only)
- Changing pysam 0-based coordinates to Bismark 1-based coordinates
- Runtime import of `template/DBiT-spatial-Methylation`

## Current state

TAPS and TAPS v2 share `scripts/call.py`, `scripts/methy_caller.py`, and
`scripts/methy_caller_CH.py`. Host calling iterates `split_bams/**/*.sorted.bam`
with `call_jobs` concurrent spot processes. Each invocation is single-threaded,
pre-enumerates reference sites, and defaults to MAPQ 1 and `max_depth` 250.
CG and CH are separate pileup passes.

The reference implementation is
`template/DBiT-spatial-Methylation/script/steps/python/06.methy_caller.py`
(host, per-spot via CB) and `06.spike_caller.py` (aggregate).

## Architecture

Stage order is unchanged:

`fastp_split -> demux_extract_bc -> align -> pool -> split -> mbias -> call -> saturation -> summary`

`split` still writes per-spot BAMs and `split_bams/per_spot_read_counts.tsv`
for saturation and summary. Host `call` no longer reads those BAMs.

```text
pooled/pooled.byCB.bam              -> methy_caller --mode host
split_barcodes (same file as split) ->   interval-parallel pileup
split_bams/per_spot_read_counts.tsv ->   CB -> XX_YY (split_bams.py rules)
call_trimming / mbias cutoffs       ->   one-pass CG and/or CH
                                        -> coverage/host/<XX>/<XX>_<YY>.CG.cov
                                        -> coverage/host/<XX>/<XX>_<YY>.CH.cov

pooled/pooled.<spike>.sorted.bam    -> methy_caller --mode aggregate
qc/mbias/host.subsampled.sorted.bam -> coverage/host_mito.*.cov
                                        coverage/<spike>.*.cov
```

### Components

| Unit | Role |
|------|------|
| `scripts/methy_caller.py` | Only caller implementation. `--mode host` reads pooled host BAM + barcodes + spot roster and writes per-spot `.cov`. `--mode aggregate` writes one `.cov` (or CG+CH pair) with no spot split. |
| `scripts/call.py` | Orchestrates host / mito / spike, forwards trim, context mode, and jobs. Host input is `pooled/pooled.byCB.bam`. |
| `scripts/methy_caller_CH.py` | Delete. |
| `scripts/make_cmd.py` and `scripts-v2/make_cmd.py` | Pass `split_barcodes` into `call.py`; drop `call_ch_caller_script`; update defaults and input checks. |

CB resolution must match `scripts/split_bams.py`, not the template spot ID
scheme:

- tag `CB` (configurable, default `CB`)
- `CB:Z:<x>+<y>`, or concatenated `xy` of length `2 * barcode_len`
- whitelist: first column, uppercase, index is 0-based order
- output stem: `{x_index:02d}_{y_index:02d}` under `coverage/host/{x_index:02d}/`

Reads with a missing tag, unparsable CB, or barcode absent from the whitelist
are skipped (same as split). Do not raise at write time.

### Parallelism

`call_jobs` means interval workers inside one chromosome (default remains 8;
the meaning changes). Implementation follows the template: `batch_size`
(default 10_000_000 bp) is reduced per chromosome so there is at least one
interval per worker; workers are a `fork` `ProcessPoolExecutor`; pileup uses
`truncate=True`.

Slurm layout is unchanged: one host sbatch and one sbatch per spike. Interval
parallelism stays inside the job. `slurm.call.host` / `slurm.call.spike`
`cpus_per_task` should be set to at least `call_jobs`.

Do not restore one sbatch per spot.

## Calling semantics (TAPS)

Host, mito, and spike use the same classification. Converted bases are
methylated (TAPS polarity). There is no EM-seq inversion in this caller.

### Shared filters

- Keep properly paired primary flags `{99, 147, 83, 163}` only.
- Skip deletion, refskip, missing `query_position`, missing sequence.
- Pileup: `ignore_overlaps=True`, `ignore_orphans=True`, `stepper="samtools"`,
  `redo_baq=False`, `compute_baq=True`.
- Missing `YD` tag: keep. `YD == "u"`: skip. No-op on typical TAPS `bwa mem`
  BAMs.
- R1/R2 left/right trim uses the existing cycle definition in
  `scripts/methy_caller.py`. CG trims both dinucleotide query positions. CH
  trims the cytosine query position; plus-strand CH also trims the neighbor
  query position.
- Defaults: `min_base_quality=30`, `min_mapping_quality=10`,
  `max_depth=50000`.
- Optional `--sample-size N` still limits to the first N reference sites per
  chromosome (debug only).

Use `query_sequence` (reference-oriented BAM SEQ). Do not use
`get_forward_sequence()` for classification.

### CpG

Pileup at the reference C of each `CG`. Observation is the dinucleotide at
`qpos` and `qpos+1`:

- `TG` or `CA` → methylated
- `CG` → unmethylated
- anything else → ignore

All four paired flags may contribute. Output line (coverage > 0 only):

```text
chrom  pos  pos  pct  meth  unmeth
```

`pos` remains pysam 0-based. `pct` is `round(100 * meth / (meth + unmeth), 2)`
with two decimal places, matching the current `format_result_line`.

### CH

Reference sites:

- Plus: reference dinucleotide in `{CA, CC, CT}`; count at the C.
- Minus: reference base `G` whose previous base maps `{T: CA, G: CC, A: CT}`;
  count at the G.

Strand flags:

- Plus: `{99, 147}` only
- Minus: `{83, 163}` only

Bases on `query_sequence`:

- Plus: `C` unmethylated, `T` methylated
- Minus: `G` unmethylated, `A` methylated

Neighbor check:

- Plus: `qpos + 1`; minus: `qpos - 1` (skip if neighbor is out of the read)
- `CA` / `CT`: neighbor must equal the expected genomic base (`A`/`T` on plus;
  `T`/`A` on minus via the reverse map)
- `CC`: plus neighbor in `{C, T}`; minus neighbor in `{G, A}`

Output line (coverage > 0 only), same contract as today:

```text
chrom  pos  pos  pct  meth  unmeth  context  strand
```

`context` is `CA`, `CC`, or `CT`. `strand` is `+` or `-`. Do not split into
separate CA/CC/CT files.

### One pileup

`call_context_mode` remains `cg` | `ch` | `both`. For `both`, one interval
pileup classifies both contexts and writes both file types. Do not run two
genome passes.

## Outputs and compatibility

### Host files

Final host files keep the current paths and columns. Interval workers write
part files that include a temporary `spot_id` column. The parent shuffles
those parts into per-spot files and drops the extra column.

Write a `.cov` for every spot listed in `split_bams/per_spot_read_counts.tsv`.
If a listed spot has no passing methylation observations, write an empty file.
Do not create files for whitelist spots that split never emitted.

Smoke mode: split emits at most 16 spots. Host calling must count and write
only those roster spots. Other CBs in the pooled BAM are ignored.

### Aggregate files

Unchanged paths:

- `coverage/host_mito.CG.cov` and optional `.CH.cov`
- `coverage/<spike_name>.CG.cov` and optional `.CH.cov`

Aggregate mode does not read barcodes or the spot roster.

Mito input remains `qc/mbias/host.subsampled.sorted.bam` when present, else
the existing subsample fallback from `pooled/pooled.byCB.bam`.

### Atomic writes and reruns

Write through a process-specific temp path and `os.replace` onto the final
path. Reruns skip work per final output file: a host roster spot whose
required `.cov` files already exist is omitted from counting and writing;
an existing mito or spike `.cov` is not recomputed. Do not skip the whole
host call because some spots are already done.

## Orchestration and config

`scripts/call.py` host mode:

- Require `--barcodes` (same TSV as `split --barcodes` / `split_barcodes`).
- Require `pooled/pooled.byCB.bam` and a `.bai` or `.csi` index.
- Require `split_bams/per_spot_read_counts.tsv`.
- Stop discovering `split_bams/**/*.sorted.bam` for host calling.

Remove `--ch-caller-script` and workflow key `call_ch_caller_script`.
`call_caller_script` still points at `scripts/methy_caller.py`.
`call_context_mode` is unchanged.

TAPS and TAPS v2 `make_cmd` must:

- Pass `--barcodes` from `split_barcodes`.
- Treat `pooled/pooled.byCB.bam`, its index, the barcode TSV, and
  `per_spot_read_counts.tsv` as host-call inputs in workdir checks (with the
  existing `--skip-workdir-input-checks` exception for `--stage all`).
- Change defaults: `call_min_mapping_quality=10`, `call_max_depth=50000`.
- Stop emitting `call_ch_caller_script`.

Keep: `call_mode`, `call_trimming` (host vs spike), `call_min_base_quality`,
`call_batch_size`, `call_chromosomes`, `call_mito_chromosomes`,
`call_sample_size`, `call_jobs`.

## Errors and dry-run

Non-zero exit when:

- BAM, BAM index, reference FASTA, or `.fai` is missing
- requested chromosomes are absent from the reference
- `--jobs`, `--batch-size`, or `--max-depth` is not `> 0`
- quality or trim values are negative
- host mode is missing `--barcodes` or `per_spot_read_counts.tsv`
- barcode whitelist is empty or has unequal lengths

`--dry-run` prints resolved parameters (mode, BAM, reference, chromosomes,
context, jobs, trims, barcodes path when host) and writes nothing.

## Documentation (same change)

Update these pages; do not copy the same default twice:

- `docs/developers/contracts.md` (`call` inputs: pooled host BAM, barcodes,
  spot roster; host no longer reads split BAMs)
- `docs/developers/config-reference.md` (new defaults, `call_jobs` meaning,
  remove `call_ch_caller_script`)
- `docs/users/taps.md` if user-facing call inputs or defaults change
- `docs/maintenance/changelog.md` under `## Current`, with `BREAKING`
- `docs/maintenance/taps-regression.md` (call `--help` / `--dry-run` still;
  drop CH-script-only checks)

Workflow JSON that currently pin `call_min_mapping_quality: 1` or
`call_max_depth: 250` should follow the new defaults unless a run explicitly
needs the old values.

`doc/` is not updated unless a task targets it.

## Verification

Required before calling the change done:

- `pixi run python scripts/make_cmd.py --version` and
  `pixi run python scripts-v2/make_cmd.py --version`
- `pixi run python scripts/methy_caller.py --help`
- `pixi run python scripts/call.py --help`
- TAPS and TAPS v2 `make_cmd --stage call --dry-run`, and `--stage all
  --dry-run` on the TAPS regression config

Unit tests (new, small synthetic sequences / fake alignments where needed):

- CG: `TG`/`CA` methylated, `CG` unmethylated, other dinucleotides ignored
- CH plus/minus flag sets
- CH CC neighbor allows C/T (plus) and G/A (minus)
- CH CA/CT neighbor must match exactly
- `YD:u` skipped; missing `YD` kept
- CB `x+y` and concatenated forms match `split_bams.py` indices and `XX_YY`
- host smoke roster: CB outside the roster is not counted
- `--dry-run` does not write outputs

A 500 Gbp benchmark is not required for this change.

## Breaking changes

Record explicitly in the changelog body:

- Host `call` input is `pooled/pooled.byCB.bam`, not per-spot split BAMs.
- `call_jobs` is interval workers, not concurrent spot processes.
- Default MAPQ 1 → 10 and `max_depth` 250 → 50000.
- CH neighbor and strand-flag rules change; CH percentages will move more
  than CG.
- CG percentages can move because of MAPQ (and any reads previously counted
  at MAPQ 1–9).
- `scripts/methy_caller_CH.py` and `call_ch_caller_script` are removed.

## Out of scope follow-ups

- Skip `split` when the user only needs calling (would require another source
  of per-spot read counts).
- Emit template-style combined `host.CG.cov` as the primary artifact.
- Port this caller to EMSeq (`--assay emseq` inversion + biscuit BAMs).
