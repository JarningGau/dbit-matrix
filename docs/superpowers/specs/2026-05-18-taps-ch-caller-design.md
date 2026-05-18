# TAPS CH Caller Design

## Summary

Add optional CH methylation calling support to the DBiT-TAPS `call` stage.

This iteration only focuses on generating `.CH.cov` outputs. It does not change
`summary`, `saturation`, `aggregate`, or any downstream consumer.

The existing CpG calling path remains the default and must keep its current
behavior and output contract.

## Goals

- Support CH methylation calling for DBiT-TAPS data.
- Keep the existing `.CG.cov` path unchanged by default.
- Allow the `call` stage to emit `.CH.cov` in addition to `.CG.cov`.
- Reuse the current caller-style output format as much as possible.
- Keep the implementation isolated so CpG logic and CH logic do not share a
  fragile code path.

## Non-Goals

- No changes to `summary`.
- No changes to `saturation`.
- No changes to `aggregate`.
- No changes to EMSeq.
- No main-pipeline validation beyond dry-run coverage for the `call` stage.

## User-Facing Behavior

The TAPS `call` stage gains a context mode switch:

- `cg`: emit only existing `.CG.cov` outputs
- `ch`: emit only `.CH.cov` outputs
- `both`: emit both `.CG.cov` and `.CH.cov` outputs

Default is `cg` so existing workflows remain unchanged.

## Configuration

Add these TAPS calling fields:

- `call_context_mode`
  - allowed values: `cg`, `ch`, `both`
  - default: `cg`
- `call_ch_caller_script`
  - default: `scripts/methy_caller_CH.py`

Existing `call_caller_script` remains the CpG caller path and continues to
default to `scripts/methy_caller.py`.

## Stage Integration

`scripts/call.py` remains the orchestration layer.

Responsibilities:

- Parse the new context mode.
- Keep current CpG command construction unchanged when `cg` is requested.
- Add parallel CH command construction when `ch` is requested.
- When mode is `both`, schedule both callers for the same BAM input.

Output naming:

- host per-spot CpG: `coverage/host/<X_index>/<X_index>_<Y_index>.CG.cov`
- host per-spot CH: `coverage/host/<X_index>/<X_index>_<Y_index>.CH.cov`
- host mito CpG: `coverage/host_mito.CG.cov`
- host mito CH: `coverage/host_mito.CH.cov`
- spike-in CpG: `coverage/<spike_name>.CG.cov`
- spike-in CH: `coverage/<spike_name>.CH.cov`

This preserves the current CpG naming contract and adds CH outputs as a parallel
artifact family.

## CH Caller Boundaries

Add a new script:

- `scripts/methy_caller_CH.py`

This script is independent from `scripts/methy_caller.py`.

Reason:

- CpG and CH have different target-site discovery rules.
- CH requires a different output schema.
- Isolating the CH path reduces regression risk for the already-stable CpG path.

## CH Output Contract

CH output is headerless and Bismark-like for the first six columns, with one
extra trailing column:

`chrom  start  end  methylation_percent  mc_count  unmeth_count  context`

Rules:

- `start` and `end` use the same coordinate convention as the existing CpG
  caller.
- `context` is one of `CA`, `CC`, or `CT`.
- Only covered sites are written.

The existing `.CG.cov` format remains unchanged.

## CH Biology And Counting Rule

For non-strand-specific DBiT-TAPS data, CH methylation is quantified by
collapsing evidence from both DNA strands onto a reference-strand CH site.

For each reference-strand CH site:

- methylated count: `NTH + DAN`
- unmethylated count: `NCH + DGN`
- total CH coverage: `NCH + NTH + DGN + DAN`
- methylation percent: `(NTH + DAN) / (NCH + NTH + DGN + DAN) * 100`

Reference:

- `feat-support-taps-CH.md`

## Target Site Discovery

`scripts/methy_caller_CH.py` scans the reference sequence and enumerates only
reference-strand CH sites:

- `CA`
- `CC`
- `CT`

Each discovered site produces at most one output row.

The row `context` is derived directly from the reference dinucleotide at that
site.

## Pileup Model

For each target CH site:

1. Run pileup over the target coordinates using the same general filtering shape
   as the current TAPS caller.
2. Apply the same read-level trimming model already used by the CpG caller:
   separate `R1` and `R2` left/right trimming.
3. Merge evidence from both orientations into the same reference-site result.

The CH caller should follow the existing TAPS caller conventions for:

- minimum base quality
- minimum mapping quality
- maximum depth
- batch size
- optional chromosome subsets
- optional dry-run mode

## Context-Specific Output Semantics

Each output row is keyed by a reference-strand cytosine in `CA`, `CC`, or `CT`.

The `context` column describes that reference-strand CH identity, not the read
orientation. Reverse-strand evidence is collapsed into the same row instead of
producing separate strand-specific records.

## Validation Scope

Required validation for this design:

- relevant CLI `--help`
- relevant workflow `--dry-run`
- functional validation of `scripts/methy_caller_CH.py` using:
  `work/test-DNAme-TAPS/pooled/pooled.puc19.sorted.bam`

This iteration does not require full host or host-mito functional validation in
the main workflow.

## Documentation Impact

If implementation changes the exposed configuration or stage contract, update:

- `docs/developers/config-reference.md`
- `docs/developers/contracts.md`

No user-guide expansion is required in this iteration because the feature is
scoped to `call` output generation only.

## Risks And Mitigations

Risk: CH counting semantics are easy to mix up with CpG semantics.
Mitigation: keep CH in a dedicated caller script with an explicit counting model.

Risk: downstream tools may assume all `.cov` files are six-column CpG outputs.
Mitigation: this iteration does not wire `.CH.cov` into downstream consumers.

Risk: introducing CH support could accidentally alter current CpG behavior.
Mitigation: keep `cg` as the default mode and preserve the existing CpG caller
path unchanged.
