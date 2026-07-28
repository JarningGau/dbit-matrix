#!/usr/bin/env python3
"""Generate and optionally submit TAPS v2 workflow commands.

TAPS v2 is a parallel driver for the same mainline stage order as
``scripts/make_cmd.py``, with methylated-linker demux:

- C→N-masked ``insert_left`` locate
- optional all-T conversion filter (``--require-c-all-t``)
- mC→T conversion QC in demux stats and sample summary

Unchanged stages reuse implementations under ``scripts/``. Demux and summary
use ``scripts-v2/extract_bc.py`` and ``scripts-v2/summary.py``.
"""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).resolve().parent.parent / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from _version import __version__
import workflow_input_checks as wic

STAGE_SEQUENCE = [
    "fastp_split",
    "demux_extract_bc",
    "align",
    "pool",
    "split",
    "mbias",
    "call",
    "saturation",
    "summary",
]
METHSCAN_STAGES = (
    "methscan_prepare",
    "methscan_filter",
    "methscan_profile",
    "methscan_smooth",
    "methscan_scan",
    "methscan_matrix",
    "methscan_all",
)
STAGE_CHOICES = [*STAGE_SEQUENCE, "aggregate", *METHSCAN_STAGES, "all"]
# Stages that may appear as top-level keys under workflow `slurm` (nested mode).
SLURM_NEST_STAGE_KEYS = frozenset(STAGE_SEQUENCE) | {"aggregate", *METHSCAN_STAGES}
STAGE_REQUIRED_FIELDS = {
    "fastp_split": ["r1", "r2", "number_of_split_parts"],
    "demux_extract_bc": [
        "barcode1_whitelist",
        "barcode2_whitelist",
        "number_of_split_parts",
    ],
    "align": ["bwa_index", "number_of_split_parts"],
    "pool": [],
    "split": ["split_barcodes"],
    "mbias": [],
    "call": ["call_reference_file", "call_chromosomes"],
    "saturation": [],
    "summary": [],
    "aggregate": [],
    "methscan_prepare": [],
    "methscan_filter": [],
    "methscan_profile": [],
    "methscan_smooth": [],
    "methscan_scan": [],
    "methscan_matrix": [],
    "methscan_all": [],
}

# (methscan_run subcommand, script basename without extension, Slurm job name prefix)
TAPS_METHSCAN_STAGE_MAP = {
    "methscan_prepare": ("prepare", "11_methscan_prepare", "dbit_methscan_prepare"),
    "methscan_filter": ("filter", "12_methscan_filter", "dbit_methscan_filter"),
    "methscan_profile": ("profile", "13_methscan_profile", "dbit_methscan_profile"),
    "methscan_smooth": ("smooth", "14_methscan_smooth", "dbit_methscan_smooth"),
    "methscan_scan": ("scan", "15_methscan_scan", "dbit_methscan_scan"),
    "methscan_matrix": ("matrix", "16_methscan_matrix", "dbit_methscan_matrix"),
    "methscan_all": ("all", "17_methscan_all", "dbit_methscan_all"),
}


def resolve_env_executable(name: str) -> str:
    candidate = Path(sys.executable).resolve().parent / name
    if candidate.is_file():
        return str(candidate)
    return name


def normalize_executable_setting(value: str | None, default_name: str) -> str:
    if not value or value == default_name:
        return resolve_env_executable(default_name)
    return value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate executable command scripts for DBiT TAPS v2 workflow "
            "(methylated-linker demux)."
        )
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    parser.add_argument(
        "--workflow-config",
        help="JSON config path for workflow/sample settings.",
    )
    parser.add_argument(
        "--runner",
        choices=["local", "slurm"],
        help="Command target: local shell or slurm sbatch.",
    )
    parser.add_argument(
        "--stage",
        choices=STAGE_CHOICES,
        help="Workflow stage to generate command script for. Default: fastp_split.",
    )
    parser.add_argument("--sample-id", help="Sample identifier.")
    parser.add_argument("--r1", help="Input R1 FASTQ(.gz).")
    parser.add_argument("--r2", help="Input R2 FASTQ(.gz).")
    parser.add_argument("--work-root", help="Work root directory. Default: work.")
    parser.add_argument(
        "--fastp-threads",
        type=int,
        help="Thread count for fastp split. Default: 8.",
    )
    parser.add_argument(
        "--number-of-split-parts",
        type=int,
        help="Value passed to fastp --split.",
    )
    parser.add_argument(
        "--fastp-bin",
        help=(
            "fastp executable path or command name. "
            "Default: fastp from current Python env if available, else fastp."
        ),
    )
    parser.add_argument(
        "--barcode1-whitelist",
        help="Whitelist path for barcodeA, used by demux stage.",
    )
    parser.add_argument(
        "--barcode2-whitelist",
        help="Whitelist path for barcodeB, used by demux stage.",
    )
    parser.add_argument(
        "--linker-bc",
        help="Linker sequence between barcode2 and barcode1 for demux stage.",
    )
    parser.add_argument(
        "--insert-left",
        help="Insert-left anchor sequence for demux stage (equivalent to previous --tn5).",
    )
    parser.add_argument("--linker1", help="Linker1 sequence for demux stage.")
    parser.add_argument("--linker2", help="Linker2 sequence for demux stage.")
    parser.add_argument("--tn5", help="Tn5 sequence for demux stage.")
    parser.add_argument(
        "--linker-edit-distance",
        type=int,
        help="Max edit distance for linker/Tn5 fallback in demux stage.",
    )
    parser.add_argument(
        "--barcode-hamming-distance",
        type=int,
        help="Max Hamming distance for whitelist fallback in demux stage.",
    )
    parser.add_argument(
        "--gzip-level",
        type=int,
        help="gzip compress level (0-9) for demux output FASTQ.",
    )
    parser.add_argument(
        "--spike-edit-distance",
        type=int,
        help=(
            "Max mismatches on non-N positions of the C→N-masked insert_left "
            "(methylated-linker) in demux. Default: 1."
        ),
    )
    parser.add_argument(
        "--require-c-all-t",
        help=(
            "Comma-separated 0-based C positions in insert_left that must all "
            "be T to keep a read (e.g. 3,6,10). Empty = no conversion filter."
        ),
    )
    parser.add_argument(
        "--bwa-index",
        help="bwa index prefix for align stage.",
    )
    parser.add_argument(
        "--bwa-threads",
        type=int,
        help="Thread count for bwa mem in align stage. Default: 2.",
    )
    parser.add_argument(
        "--bwa-bin",
        help=(
            "bwa executable path or command name. "
            "Default: bwa from current Python env if available, else bwa."
        ),
    )
    parser.add_argument(
        "--sinto-bin",
        help=(
            "sinto executable path or command name. "
            "Default: sinto from current Python env if available, else sinto."
        ),
    )
    parser.add_argument(
        "--samtools-bin",
        help=(
            "samtools executable path or command name. "
            "Default: samtools from current Python env if available, else samtools."
        ),
    )
    parser.add_argument(
        "--samtools-threads",
        type=int,
        help="Thread count for samtools sort in pool stage. Default: 4.",
    )
    parser.add_argument(
        "--host-sort-mem",
        help="Memory per thread for host samtools sort -m in pool stage. Default: 16G.",
    )
    parser.add_argument(
        "--split-barcodes",
        help="Barcode whitelist path for split stage. Default: barcode1_whitelist.",
    )
    parser.add_argument(
        "--split-cb-tag",
        help="CB-like tag name used in split stage. Default: CB.",
    )
    parser.add_argument(
        "--split-threads-read",
        type=int,
        help="Reader threads for split stage. Default: 1.",
    )
    parser.add_argument(
        "--split-threads-write",
        type=int,
        help="Writer threads for split stage. Default: 1.",
    )
    parser.add_argument(
        "--split-smoke",
        action="store_true",
        default=None,
        help="Enable smoke mode for split stage and emit up to 16 spot BAMs.",
    )
    parser.add_argument(
        "--split-sort-jobs",
        type=int,
        help="Parallel job count for bam_sort_parallel in split stage. Default: 8.",
    )
    parser.add_argument(
        "--mbias-host-subsample-fraction",
        type=float,
        help="Host subsampling fraction for mbias stage. Default: 0.1.",
    )
    parser.add_argument(
        "--mbias-mode",
        choices=["all", "host", "spike"],
        help="Mode for mbias stage. Default: spike.",
    )
    parser.add_argument(
        "--mbias-max-cycle",
        type=int,
        help="Maximum read cycle for mbias stage. Default: 150.",
    )
    parser.add_argument(
        "--mbias-min-mapping-quality",
        type=int,
        help="Minimum mapping quality for mbias stage. Default: 1.",
    )
    parser.add_argument(
        "--mbias-script",
        help="Path to mbias script. Default: scripts/mbias.py.",
    )
    parser.add_argument(
        "--call-reference-file",
        help="Reference FASTA path for call stage.",
    )
    parser.add_argument(
        "--call-chromosomes",
        help="Comma-separated chromosomes for call stage output (e.g. chr1,chr2).",
    )
    parser.add_argument(
        "--call-mito-chromosomes",
        help="Comma-separated chromosomes for host_mito outputs. Default: chrM.",
    )
    parser.add_argument(
        "--call-jobs",
        type=int,
        help="Parallel host spot jobs for call stage. Default: 8.",
    )
    parser.add_argument(
        "--call-min-base-quality",
        type=int,
        help="min_base_quality passed to methy_caller in call stage. Default: 30.",
    )
    parser.add_argument(
        "--call-min-mapping-quality",
        type=int,
        help="min_mapping_quality passed to methy_caller in call stage. Default: 1.",
    )
    parser.add_argument(
        "--call-sample-size",
        type=int,
        help="Optional sample_size passed to methy_caller in call stage.",
    )
    parser.add_argument(
        "--call-max-depth",
        type=int,
        help="max_depth passed to methy_caller in call stage. Default: 250.",
    )
    parser.add_argument(
        "--call-batch-size",
        type=int,
        help="batch_size passed to methy_caller in call stage. Default: 10000000.",
    )
    parser.add_argument(
        "--call-mode",
        choices=["all", "host", "spike"],
        help="Mode for call stage. Default: all.",
    )
    parser.add_argument(
        "--call-context-mode",
        choices=["cg", "ch", "both"],
        help="Context mode for TAPS call stage. Default: cg.",
    )
    parser.add_argument(
        "--call-r1-left-trimming",
        type=int,
        help="R1 left-end trimming for call stage. Default: 0.",
    )
    parser.add_argument(
        "--call-r1-right-trimming",
        type=int,
        help="R1 right-end trimming for call stage. Default: 0.",
    )
    parser.add_argument(
        "--call-r2-left-trimming",
        type=int,
        help="R2 left-end trimming for call stage. Default: 0.",
    )
    parser.add_argument(
        "--call-r2-right-trimming",
        type=int,
        help="R2 right-end trimming for call stage. Default: 0.",
    )
    parser.add_argument(
        "--call-caller-script",
        help="Path to methy_caller script for call stage. Default: scripts/methy_caller.py.",
    )
    parser.add_argument(
        "--call-ch-caller-script",
        help="Path to CH caller script for call stage. Default: scripts/methy_caller_CH.py.",
    )
    parser.add_argument(
        "--saturation-script",
        help="Path to saturation script. Default: scripts/saturation.py.",
    )
    parser.add_argument(
        "--saturation-reads-threshold",
        type=float,
        help="HQ spot reads threshold for saturation stage. Default: 1e6.",
    )
    parser.add_argument(
        "--summary-script",
        help="Path to summary script. Default: scripts-v2/summary.py.",
    )
    parser.add_argument(
        "--aggregate-script",
        help="Path to aggregate script. Default: scripts/aggregate.py.",
    )
    parser.add_argument(
        "--aggregate-sort-mem",
        help=(
            "Memory limit for GNU `sort -S` in aggregate stage "
            "(passed to scripts/aggregate.py via `--sort-mem`). "
            "Default comes from workflow config `aggregate_sort_mem` or 8G."
        ),
    )
    parser.add_argument(
        "--methscan-pixi-manifest",
        help=(
            "Directory with pixi.toml for methscan (passed to --pixi-manifest). "
            "Default: envs/methscan under repo root."
        ),
    )
    parser.add_argument(
        "--methscan-run-script",
        help="Path to methscan_run.py. Default: scripts/methscan_run.py.",
    )
    parser.add_argument(
        "--methscan-prepare-chunksize",
        type=int,
        help="Passed to methscan prepare --chunksize (default: 10000000).",
    )
    parser.add_argument(
        "--methscan-filter-min-sites",
        type=int,
        help="Passed to methscan filter --min-sites (default: 50000).",
    )
    parser.add_argument(
        "--methscan-tss-bed",
        help="BED file for methscan profile / methscan_all (required for those stages).",
    )
    parser.add_argument(
        "--methscan-profile-strand-column",
        type=int,
        help="methscan profile --strand-column (default: 6).",
    )
    parser.add_argument(
        "--methscan-profile-prepared-subdir",
        choices=("compact", "filter"),
        help="methscan profile --prepared-dir (default: compact).",
    )
    parser.add_argument(
        "--methscan-profile-csv",
        help="methscan profile output CSV path (optional).",
    )
    parser.add_argument(
        "--methscan-smooth-bandwidth",
        type=int,
        help="methscan smooth --bandwidth (optional).",
    )
    parser.add_argument(
        "--methscan-smooth-use-weights",
        action="store_true",
        default=None,
        help=(
            "Pass --use-weights to methscan smooth. "
            "Omit to take methscan_smooth_use_weights from workflow JSON."
        ),
    )
    parser.add_argument(
        "--methscan-scan-threads",
        type=int,
        help="methscan scan --threads (default: 10).",
    )
    parser.add_argument(
        "--methscan-vmrs-bed",
        help="methscan scan/matrix VMRs .bed path (optional; default methscan/VMRs.bed).",
    )
    parser.add_argument(
        "--methscan-matrix-threads",
        type=int,
        help="methscan matrix --threads (default: 10).",
    )
    parser.add_argument(
        "--methscan-matrix-sparse",
        action="store_true",
        default=None,
        help=(
            "Pass --sparse to methscan matrix. "
            "Omit to take methscan_matrix_sparse from workflow JSON."
        ),
    )
    parser.add_argument(
        "--methscan-matrix-prefix",
        help="methscan matrix output directory (optional).",
    )
    parser.add_argument(
        "--spike-in-index",
        action="append",
        default=None,
        help=(
            "Spike-in reference in NAME=INDEX format. "
            "May be specified multiple times."
        ),
    )
    parser.add_argument(
        "--submit",
        action="store_true",
        help="Submit immediately after generating command file.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print command and output path without writing files.",
    )
    parser.add_argument(
        "--skip-workdir-input-checks",
        action="store_true",
        help=(
            "Do not require prior-stage outputs under the sample work directory "
            "(shard_fastq, demux shards, align_shards, pooled BAM). "
            "When generating --stage all, this is passed to each per-stage subprocess automatically."
        ),
    )

    # Slurm-only options.
    parser.add_argument("--slurm-partition")
    parser.add_argument("--slurm-mem")
    parser.add_argument("--slurm-cpus-per-task", type=int)
    parser.add_argument("--slurm-output")
    parser.add_argument("--slurm-error")
    return parser.parse_args()


def quoted(parts: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in parts)


def build_fastp_split_command(args: argparse.Namespace, sample_work: Path) -> str:
    script_path = Path("scripts/fastp_split.py")
    command = [
        sys.executable,
        str(script_path),
        "--r1",
        args.r1,
        "--r2",
        args.r2,
        "--work-path",
        str(sample_work),
        "--fastp-threads",
        str(args.fastp_threads),
        "--number-of-split-parts",
        str(args.number_of_split_parts),
        "--fastp-bin",
        args.fastp_bin,
    ]
    return quoted(command)


def build_demux_chunk_command(
    args: argparse.Namespace, r1_path: Path, r2_path: Path, out_prefix: Path
) -> str:
    script_path = Path("scripts-v2/extract_bc.py")
    return quoted(
        [
            sys.executable,
            str(script_path),
            str(r1_path),
            str(r2_path),
            "-b1",
            args.barcode1_whitelist,
            "-b2",
            args.barcode2_whitelist,
            "-o",
            str(out_prefix),
            "--linker-bc",
            args.linker_bc,
            "--insert-left",
            args.insert_left,
            "--linker-edit-distance",
            str(args.linker_edit_distance),
            "--barcode-hamming-distance",
            str(args.barcode_hamming_distance),
            "--gzip-level",
            str(args.gzip_level),
            "--spike-edit-distance",
            str(args.spike_edit_distance),
            "--require-c-all-t",
            args.require_c_all_t,
        ]
    )


def build_demux_local_batch_command(args: argparse.Namespace, sample_work: Path) -> str:
    chunk_dir = sample_work / "shard_fastq"
    demux_dir = sample_work / "demux"
    chunk_dir_q = shlex.quote(str(chunk_dir))
    demux_dir_q = shlex.quote(str(demux_dir))
    py = quoted([sys.executable, "scripts-v2/extract_bc.py"])
    return (
        f"chunk_dir={chunk_dir_q}\n"
        f"demux_dir={demux_dir_q}\n"
        "bar_width=30\n"
        "\n"
        "mkdir -p \"$demux_dir\"\n"
        "shopt -s nullglob\n"
        "r1_files=(\"$chunk_dir\"/*.R1.fq.gz)\n"
        "total=${#r1_files[@]}\n"
        "if [ \"$total\" -eq 0 ]; then\n"
        '  echo "[demux] no chunk found under shard_fastq"\n'
        "  exit 1\n"
        "fi\n"
        "\n"
        "idx=0\n"
        "for r1 in \"${r1_files[@]}\"; do\n"
        "  idx=$((idx + 1))\n"
        "  percent=$((idx * 100 / total))\n"
        "  filled=$((idx * bar_width / total))\n"
        "  empty=$((bar_width - filled))\n"
        "  bar_filled=$(printf '%*s' \"$filled\" '' | tr ' ' '#')\n"
        "  bar_empty=$(printf '%*s' \"$empty\" '')\n"
        '  chunk="$(basename "$r1" .R1.fq.gz)"\n'
        "  r2=\"$chunk_dir/${chunk}.R2.fq.gz\"\n"
        "  printf '[demux] [%s%s] %3d%% (%d/%d) %s\\n' "
        "\"$bar_filled\" \"$bar_empty\" \"$percent\" \"$idx\" \"$total\" \"$chunk\"\n"
        '  [ -f "$r2" ] || { echo "[demux] missing pair for $r1: $r2"; exit 1; }\n'
        f"  {py} "
        '"$r1" "$r2" '
        f"-b1 {shlex.quote(args.barcode1_whitelist)} "
        f"-b2 {shlex.quote(args.barcode2_whitelist)} "
        f"-o \"$demux_dir\"/${{chunk}} "
        f"--linker-bc {shlex.quote(args.linker_bc)} "
        f"--insert-left {shlex.quote(args.insert_left)} "
        f"--linker-edit-distance {int(args.linker_edit_distance)} "
        f"--barcode-hamming-distance {int(args.barcode_hamming_distance)} "
        f"--gzip-level {int(args.gzip_level)} "
        f"--spike-edit-distance {int(args.spike_edit_distance)} "
        f"--require-c-all-t {shlex.quote(args.require_c_all_t)}\n"
        "done\n"
        "\n"
        'echo "[demux] done"'
    )


def build_align_chunk_command(
    args: argparse.Namespace, sample_work: Path, chunk: str
) -> str:
    script_path = Path("scripts/align.py")
    command = [
        sys.executable,
        str(script_path),
        "--work-path",
        str(sample_work),
        "--chunk",
        chunk,
        "--bwa-index",
        args.bwa_index,
        "--bwa-threads",
        str(args.bwa_threads),
        "--bwa-bin",
        args.bwa_bin,
        "--sinto-bin",
        args.sinto_bin,
        "--samtools-bin",
        args.samtools_bin,
    ]
    for item in args.spike_in_index:
        command.extend(["--spike-in-index", item])
    return quoted(command)


def build_align_local_batch_command(args: argparse.Namespace, sample_work: Path) -> str:
    script_path = Path("scripts/align.py")
    command = [
        sys.executable,
        str(script_path),
        "--work-path",
        str(sample_work),
        "--bwa-index",
        args.bwa_index,
        "--bwa-threads",
        str(args.bwa_threads),
        "--bwa-bin",
        args.bwa_bin,
        "--sinto-bin",
        args.sinto_bin,
        "--samtools-bin",
        args.samtools_bin,
    ]
    for item in args.spike_in_index:
        command.extend(["--spike-in-index", item])
    return quoted(command)


def parse_spike_names(spike_in_index: list[str]) -> list[str]:
    names: list[str] = []
    for item in spike_in_index:
        name, _, _ = item.partition("=")
        name = name.strip()
        if name and name not in names:
            names.append(name)
    return names


def build_pool_command(
    args: argparse.Namespace, sample_work: Path, mode: str
) -> str:
    script_path = Path("scripts/pool.py")
    command = [
        sys.executable,
        str(script_path),
        "--work-path",
        str(sample_work),
        "--samtools-bin",
        args.samtools_bin,
        "--samtools-threads",
        str(args.samtools_threads),
        "--host-sort-mem",
        args.host_sort_mem,
        "--mode",
        mode,
    ]
    for spike_name in parse_spike_names(args.spike_in_index):
        command.extend(["--spike-in-name", spike_name])
    return quoted(command)


def build_split_command(args: argparse.Namespace, sample_work: Path) -> str:
    script_path = Path("scripts/split_bams.py")
    command = [
        sys.executable,
        str(script_path),
        "--in-bam",
        str(sample_work / "pooled" / "pooled.byCB.bam"),
        "--barcodes",
        args.split_barcodes,
        "--out-dir",
        str(sample_work / "split_bams"),
        "--cb-tag",
        args.split_cb_tag,
        "--threads-read",
        str(args.split_threads_read),
        "--threads-write",
        str(args.split_threads_write),
    ]
    if args.split_smoke:
        command.append("--smoke")
    return quoted(command)


def build_split_sort_command(args: argparse.Namespace, sample_work: Path) -> str:
    script_path = Path("scripts/bam_sort_parallel.py")
    command = [
        sys.executable,
        str(script_path),
        "--work-path",
        str(sample_work),
        "--samtools-bin",
        args.samtools_bin,
        "--jobs",
        str(args.split_sort_jobs),
    ]
    return quoted(command)


CALL_TRIM_FIELDS = ("r1_left", "r1_right", "r2_left", "r2_right")


def resolve_call_trim_target(
    cfg_target: object,
    flat_defaults: dict[str, int],
    target_label: str,
) -> dict[str, int]:
    if cfg_target is None:
        cfg_target = {}
    if not isinstance(cfg_target, dict):
        raise ValueError(f"call_trimming.{target_label} must be an object")
    result: dict[str, int] = {}
    for field in CALL_TRIM_FIELDS:
        raw = cfg_target.get(field, flat_defaults[field])
        if raw is None:
            raw = 0
        value = int(raw)
        if value < 0:
            raise ValueError(f"call_trimming.{target_label}.{field} must be >= 0")
        result[field] = value
    return result


def build_call_command(
    args: argparse.Namespace,
    sample_work: Path,
    mode: str,
    spike_name: str | None = None,
) -> str:
    script_path = Path("scripts/call.py")
    trim = args.call_trim_spike if mode == "spike" else args.call_trim_host
    command = [
        sys.executable,
        str(script_path),
        "--work-path",
        str(sample_work),
        "--mode",
        mode,
        "--reference-file",
        args.call_reference_file,
        "--chromosomes",
        args.call_chromosomes,
        "--mito-chromosomes",
        args.call_mito_chromosomes,
        "--context-mode",
        args.call_context_mode,
        "--jobs",
        str(args.call_jobs),
        "--min-base-quality",
        str(args.call_min_base_quality),
        "--min-mapping-quality",
        str(args.call_min_mapping_quality),
        "--max-depth",
        str(args.call_max_depth),
        "--batch-size",
        str(args.call_batch_size),
        "--r1-left-trimming",
        str(trim["r1_left"]),
        "--r1-right-trimming",
        str(trim["r1_right"]),
        "--r2-left-trimming",
        str(trim["r2_left"]),
        "--r2-right-trimming",
        str(trim["r2_right"]),
        "--caller-script",
        args.call_caller_script,
        "--ch-caller-script",
        args.call_ch_caller_script,
        "--samtools-bin",
        args.samtools_bin,
        "--samtools-threads",
        str(args.samtools_threads),
        "--host-subsample-fraction",
        str(args.mbias_host_subsample_fraction),
    ]
    if args.call_sample_size is not None:
        command.extend(["--sample-size", str(args.call_sample_size)])
    if mode in ("all", "spike"):
        for item in args.spike_in_index:
            command.extend(["--spike-reference", item])
    if spike_name:
        command.extend(["--spike-in-name", spike_name])
    return quoted(command)


def build_mbias_command(
    args: argparse.Namespace,
    sample_work: Path,
    mode: str,
    spike_name: str | None = None,
) -> str:
    script_path = Path(args.mbias_script)
    command = [
        sys.executable,
        str(script_path),
        "--work-path",
        str(sample_work),
        "--mode",
        mode,
        "--samtools-bin",
        args.samtools_bin,
        "--samtools-threads",
        str(args.samtools_threads),
        "--host-subsample-fraction",
        str(args.mbias_host_subsample_fraction),
        "--max-cycle",
        str(args.mbias_max_cycle),
        "--min-mapping-quality",
        str(args.mbias_min_mapping_quality),
    ]
    if mode in ("all", "host"):
        command.extend(["--reference-file", args.call_reference_file])
        command.extend(["--chromosomes", args.call_chromosomes])
    if mode in ("all", "spike"):
        for item in args.spike_in_index:
            command.extend(["--spike-reference", item])
    if spike_name:
        command.extend(["--spike-in-name", spike_name])
    return quoted(command)


def build_aggregate_command(args: argparse.Namespace, sample_work: Path) -> str:
    script_path = Path(args.aggregate_script)
    command = [
        sys.executable,
        str(script_path),
        "--work-path",
        str(sample_work),
    ]
    sort_mem = getattr(args, "sort_mem", None)
    if sort_mem:
        command.extend(["--sort-mem", str(sort_mem)])
    return quoted(command)


def build_methscan_run_command(
    step: str,
    args: argparse.Namespace,
    sample_work: Path,
) -> str:
    script_path = Path(args.methscan_run_script)
    command = [sys.executable, str(script_path), step, "--work-path", str(sample_work)]
    manifest = getattr(args, "methscan_pixi_manifest", None)
    if manifest:
        command.extend(["--pixi-manifest", str(manifest)])
    if step == "prepare":
        command.extend(
            ["--chunksize", str(getattr(args, "methscan_prepare_chunksize", 10_000_000))]
        )
    elif step == "filter":
        command.extend(
            ["--min-sites", str(getattr(args, "methscan_filter_min_sites", 50_000))]
        )
    elif step == "profile":
        command.extend(
            [
                "--tss-bed",
                str(args.methscan_tss_bed),
                "--strand-column",
                str(getattr(args, "methscan_profile_strand_column", 6)),
                "--prepared-dir",
                str(getattr(args, "methscan_profile_prepared_subdir", "compact")),
            ]
        )
        if getattr(args, "methscan_profile_csv", None):
            command.extend(["--profile-csv", str(args.methscan_profile_csv)])
    elif step == "smooth":
        bw = getattr(args, "methscan_smooth_bandwidth", None)
        if bw is not None:
            command.extend(["--bandwidth", str(bw)])
        if getattr(args, "methscan_smooth_use_weights", False):
            command.append("--use-weights")
    elif step == "scan":
        command.extend(
            ["--threads", str(getattr(args, "methscan_scan_threads", 10))]
        )
        if getattr(args, "methscan_vmrs_bed", None):
            command.extend(["--vmrs-bed", str(args.methscan_vmrs_bed)])
    elif step == "matrix":
        command.extend(
            ["--threads", str(getattr(args, "methscan_matrix_threads", 10))]
        )
        if getattr(args, "methscan_matrix_sparse", False):
            command.append("--sparse")
        if getattr(args, "methscan_vmrs_bed", None):
            command.extend(["--vmrs-bed", str(args.methscan_vmrs_bed)])
        if getattr(args, "methscan_matrix_prefix", None):
            command.extend(["--matrix-prefix", str(args.methscan_matrix_prefix)])
    elif step == "all":
        command.extend(
            [
                "--chunksize",
                str(getattr(args, "methscan_prepare_chunksize", 10_000_000)),
                "--min-sites",
                str(getattr(args, "methscan_filter_min_sites", 50_000)),
                "--tss-bed",
                str(args.methscan_tss_bed),
                "--strand-column",
                str(getattr(args, "methscan_profile_strand_column", 6)),
                "--prepared-dir",
                str(getattr(args, "methscan_profile_prepared_subdir", "compact")),
            ]
        )
        if getattr(args, "methscan_profile_csv", None):
            command.extend(["--profile-csv", str(args.methscan_profile_csv)])
        bw = getattr(args, "methscan_smooth_bandwidth", None)
        if bw is not None:
            command.extend(["--bandwidth", str(bw)])
        if getattr(args, "methscan_smooth_use_weights", False):
            command.append("--use-weights")
        command.extend(
            [
                "--scan-threads",
                str(getattr(args, "methscan_scan_threads", 10)),
                "--matrix-threads",
                str(getattr(args, "methscan_matrix_threads", 10)),
            ]
        )
        if getattr(args, "methscan_matrix_sparse", False):
            command.append("--sparse")
        if getattr(args, "methscan_vmrs_bed", None):
            command.extend(["--vmrs-bed", str(args.methscan_vmrs_bed)])
        if getattr(args, "methscan_matrix_prefix", None):
            command.extend(["--matrix-prefix", str(args.methscan_matrix_prefix)])
    else:
        raise ValueError(f"unknown methscan step: {step}")
    return quoted(command)


def methscan_namespace_from_settings(settings: dict) -> argparse.Namespace:
    return argparse.Namespace(
        methscan_run_script=settings["methscan_run_script"],
        methscan_pixi_manifest=settings.get("methscan_pixi_manifest"),
        methscan_prepare_chunksize=settings["methscan_prepare_chunksize"],
        methscan_filter_min_sites=settings["methscan_filter_min_sites"],
        methscan_tss_bed=settings.get("methscan_tss_bed"),
        methscan_profile_strand_column=settings["methscan_profile_strand_column"],
        methscan_profile_prepared_subdir=settings["methscan_profile_prepared_subdir"],
        methscan_profile_csv=settings.get("methscan_profile_csv"),
        methscan_smooth_bandwidth=settings.get("methscan_smooth_bandwidth"),
        methscan_smooth_use_weights=settings["methscan_smooth_use_weights"],
        methscan_scan_threads=settings["methscan_scan_threads"],
        methscan_vmrs_bed=settings.get("methscan_vmrs_bed"),
        methscan_matrix_threads=settings["methscan_matrix_threads"],
        methscan_matrix_sparse=settings["methscan_matrix_sparse"],
        methscan_matrix_prefix=settings.get("methscan_matrix_prefix"),
    )


def emit_methscan_stage_command(
    *,
    settings: dict,
    sample_work: Path,
    command_dir: Path,
    log_dir: Path,
    step: str,
    script_filename: str,
    job_name: str,
    generated_scripts: list[Path],
) -> None:
    command_args = methscan_namespace_from_settings(settings)
    command = build_methscan_run_command(step, command_args, sample_work)
    script_path = command_dir / script_filename
    print(f"[make_cmd] runner={settings['runner']}")
    print(f"[make_cmd] stage={settings['stage']}")
    print(f"[make_cmd] sample_id={settings['sample_id']}")
    print(f"[make_cmd] script={script_path}")
    print(f"[make_cmd] command={command}")
    if settings["dry_run"]:
        return
    if settings["runner"] == "local":
        generate_local_script(command, script_path)
    else:
        slurm_args = argparse.Namespace(
            job_name=job_name,
            slurm_partition=settings["slurm_partition"],
            slurm_mem=settings["slurm_mem"],
            slurm_cpus_per_task=settings["slurm_cpus_per_task"],
            slurm_output=settings["slurm_output"].replace("%x", job_name),
            slurm_error=settings["slurm_error"].replace("%x", job_name),
            module_line="",
        )
        generate_slurm_script(command, script_path, log_dir, slurm_args)
    generated_scripts.append(script_path)


def build_summary_command(args: argparse.Namespace, sample_work: Path) -> str:
    script_path = Path(args.summary_script)
    command = [
        sys.executable,
        str(script_path),
        "--work-path",
        str(sample_work),
        "--mito-chromosomes",
        args.call_mito_chromosomes,
    ]
    for spike_name in parse_spike_names(args.spike_in_index):
        command.extend(["--spike-in-name", spike_name])
    return quoted(command)


def build_saturation_command(args: argparse.Namespace, sample_work: Path) -> str:
    script_path = Path(args.saturation_script)
    command = [
        sys.executable,
        str(script_path),
        "--work-path",
        str(sample_work),
        "--reads-threshold",
        str(args.saturation_reads_threshold),
    ]
    return quoted(command)


def build_chunk_names(number_of_split_parts: int) -> list[str]:
    if number_of_split_parts <= 0:
        raise ValueError("number_of_split_parts must be > 0")
    width = max(4, len(str(number_of_split_parts)))
    return [f"{index:0{width}d}" for index in range(1, number_of_split_parts + 1)]


def build_demux_chunks_from_config(
    sample_work: Path, number_of_split_parts: int
) -> list[tuple[str, Path, Path, Path]]:
    chunk_dir = sample_work / "shard_fastq"
    demux_dir = sample_work / "demux"
    return [
        (
            chunk,
            chunk_dir / f"{chunk}.R1.fq.gz",
            chunk_dir / f"{chunk}.R2.fq.gz",
            demux_dir / chunk,
        )
        for chunk in build_chunk_names(number_of_split_parts)
    ]


def build_align_chunks_from_config(
    sample_work: Path, number_of_split_parts: int
) -> list[tuple[str, Path, Path]]:
    demux_dir = sample_work / "demux"
    return [
        (
            chunk,
            demux_dir / f"{chunk}.R1.demux.fq.gz",
            demux_dir / f"{chunk}.R2.demux.fq.gz",
        )
        for chunk in build_chunk_names(number_of_split_parts)
    ]


def normalize_spike_in_index(raw) -> list[str]:
    if raw in (None, ""):
        return []
    if isinstance(raw, list):
        values = raw
    elif isinstance(raw, dict):
        values = [f"{name}={index}" for name, index in raw.items()]
    else:
        raise ValueError("spike_in_index must be a list or object in workflow config")
    normalized: list[str] = []
    for item in values:
        if not isinstance(item, str):
            raise ValueError("spike_in_index entries must be strings")
        name, sep, index = item.partition("=")
        if not sep or not name.strip() or not index.strip():
            raise ValueError(
                f"invalid spike_in_index entry (expected NAME=INDEX): {item}"
            )
        normalized.append(f"{name.strip()}={index.strip()}")
    return normalized


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def load_workflow_config(path: str) -> dict:
    config_path = Path(path)
    with config_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError("workflow config must be a JSON object")
    return data


def pick(cli_value, cfg_value):
    return cli_value if cli_value is not None else cfg_value


def validate_required_for_stage(stage: str, settings: dict) -> None:
    required = ["runner", "sample_id", *STAGE_REQUIRED_FIELDS[stage]]
    missing = [key for key in required if settings.get(key) in (None, "")]
    if missing:
        raise ValueError(f"missing required settings: {', '.join(missing)}")
    if stage in ("methscan_profile", "methscan_all") and not settings.get(
        "methscan_tss_bed"
    ):
        raise ValueError("methscan_tss_bed is required for this stage")


def validate_inputs_for_stage(
    stage: str,
    settings: dict,
    sample_work: Path,
    *,
    skip_workdir_inputs: bool = False,
) -> None:
    """Raise ValueError if required input paths are missing before script generation."""
    spike = settings["spike_in_index"]

    if stage == "fastp_split":
        wic.require_file("r1", wic.resolve_config_path(settings["r1"]))
        wic.require_file("r2", wic.resolve_config_path(settings["r2"]))
        wic.require_optional_executable_path("fastp_bin", settings["fastp_bin"])
    elif stage == "demux_extract_bc":
        wic.require_file(
            "barcode1_whitelist",
            wic.resolve_config_path(settings["barcode1_whitelist"]),
        )
        wic.require_file(
            "barcode2_whitelist",
            wic.resolve_config_path(settings["barcode2_whitelist"]),
        )
        wic.require_script_path(
            "scripts-v2/extract_bc.py",
            "scripts-v2/extract_bc.py",
        )
        if not skip_workdir_inputs:
            n = int(settings["number_of_split_parts"])
            chunk_dir = sample_work / "shard_fastq"
            for ch in wic.chunk_names(n):
                wic.require_file(
                    f"shard_fastq/{ch}.R1.fq.gz",
                    chunk_dir / f"{ch}.R1.fq.gz",
                )
                wic.require_file(
                    f"shard_fastq/{ch}.R2.fq.gz",
                    chunk_dir / f"{ch}.R2.fq.gz",
                )
    elif stage == "align":
        wic.require_bwa_index_prefix(settings["bwa_index"])
        wic.require_spike_index_paths(spike)
        if not skip_workdir_inputs:
            n = int(settings["number_of_split_parts"])
            demux_dir = sample_work / "demux"
            for ch in wic.chunk_names(n):
                wic.require_file(
                    f"demux/{ch}.R1.demux.fq.gz",
                    demux_dir / f"{ch}.R1.demux.fq.gz",
                )
                wic.require_file(
                    f"demux/{ch}.R2.demux.fq.gz",
                    demux_dir / f"{ch}.R2.demux.fq.gz",
                )
        wic.require_optional_executable_path("bwa_bin", settings["bwa_bin"])
        wic.require_optional_executable_path("sinto_bin", settings["sinto_bin"])
        wic.require_optional_executable_path("samtools_bin", settings["samtools_bin"])
    elif stage == "pool":
        if not skip_workdir_inputs:
            wic.validate_pool_align_shards(sample_work, spike)
        wic.require_optional_executable_path("samtools_bin", settings["samtools_bin"])
    elif stage == "split":
        if not skip_workdir_inputs:
            wic.require_file(
                "pooled.byCB.bam",
                sample_work / "pooled" / "pooled.byCB.bam",
            )
        wic.require_file(
            "split_barcodes",
            wic.resolve_config_path(settings["split_barcodes"]),
        )
        wic.require_optional_executable_path("samtools_bin", settings["samtools_bin"])
    elif stage == "mbias":
        wic.require_script_path("mbias_script", settings["mbias_script"])
        mbias_mode = settings.get("mbias_mode") or "spike"
        if mbias_mode in ("all", "host"):
            wic.require_file(
                "call_reference_file",
                wic.resolve_config_path(settings["call_reference_file"]),
            )
        if mbias_mode in ("all", "spike") and parse_spike_names(spike):
            wic.require_spike_index_paths(spike)
        wic.require_optional_executable_path("samtools_bin", settings["samtools_bin"])
    elif stage == "call":
        wic.require_script_path("call_caller_script", settings["call_caller_script"])
        call_mode = settings.get("call_mode") or "all"
        if call_mode in ("all", "host"):
            wic.require_file(
                "call_reference_file",
                wic.resolve_config_path(settings["call_reference_file"]),
            )
        if call_mode in ("all", "spike") and parse_spike_names(spike):
            wic.require_spike_index_paths(spike)
        wic.require_optional_executable_path("samtools_bin", settings["samtools_bin"])
    elif stage == "saturation":
        wic.require_script_path("saturation_script", settings["saturation_script"])
    elif stage == "summary":
        wic.require_script_path("summary_script", settings["summary_script"])
    elif stage == "aggregate":
        wic.require_script_path("aggregate_script", settings["aggregate_script"])
    elif stage in METHSCAN_STAGES:
        wic.require_script_path("methscan_run_script", settings["methscan_run_script"])
        wic.validate_optional_path_dir(
            "methscan_pixi_manifest", settings.get("methscan_pixi_manifest")
        )
        if stage in ("methscan_profile", "methscan_all"):
            if not settings.get("methscan_tss_bed"):
                raise ValueError("methscan_tss_bed is required for this stage")
            wic.require_file(
                "methscan_tss_bed",
                wic.resolve_config_path(settings["methscan_tss_bed"]),
            )
        wic.validate_optional_path_file(
            "methscan_profile_csv", settings.get("methscan_profile_csv")
        )
        if stage in ("methscan_scan", "methscan_matrix", "methscan_all"):
            wic.validate_optional_path_file(
                "methscan_vmrs_bed", settings.get("methscan_vmrs_bed")
            )
    else:
        raise ValueError(f"unsupported stage for input validation: {stage}")


def select_stage_slurm_cfg(slurm_cfg_raw: dict, stage: str) -> dict:
    if any(key in slurm_cfg_raw for key in SLURM_NEST_STAGE_KEYS):
        stage_slurm_cfg = slurm_cfg_raw.get(stage, {})
    else:
        # Backward compatibility: legacy flat slurm config.
        stage_slurm_cfg = slurm_cfg_raw
    if stage_slurm_cfg is None:
        stage_slurm_cfg = {}
    if not isinstance(stage_slurm_cfg, dict):
        raise ValueError("selected slurm config must be an object")
    return stage_slurm_cfg


def apply_stage_slurm_settings(
    settings: dict, args: argparse.Namespace, stage: str
) -> dict:
    stage_settings = dict(settings)
    stage_settings["stage"] = stage
    stage_slurm_cfg = select_stage_slurm_cfg(settings["slurm_cfg_raw"], stage)
    split_bams_slurm_cfg = resolve_step_slurm_cfg(
        stage_slurm_cfg,
        "split_bams",
        {"split_bams", "sort"},
    )
    split_sort_slurm_cfg = resolve_step_slurm_cfg(
        stage_slurm_cfg,
        "sort",
        {"split_bams", "sort"},
    )
    call_host_slurm_cfg = resolve_step_slurm_cfg(
        stage_slurm_cfg,
        "host",
        {"host", "spike"},
    )
    call_spike_slurm_cfg = resolve_step_slurm_cfg(
        stage_slurm_cfg,
        "spike",
        {"host", "spike"},
    )
    mbias_host_slurm_cfg = resolve_step_slurm_cfg(
        stage_slurm_cfg,
        "host",
        {"host", "spike"},
    )
    mbias_spike_slurm_cfg = resolve_step_slurm_cfg(
        stage_slurm_cfg,
        "spike",
        {"host", "spike"},
    )

    stage_settings["slurm_partition"] = (
        pick(args.slurm_partition, stage_slurm_cfg.get("partition"))
        or settings["slurm_partition"]
    )
    stage_settings["slurm_mem"] = (
        pick(args.slurm_mem, stage_slurm_cfg.get("mem")) or settings["slurm_mem"]
    )
    stage_settings["slurm_cpus_per_task"] = (
        pick(args.slurm_cpus_per_task, stage_slurm_cfg.get("cpus_per_task"))
        or settings["slurm_cpus_per_task"]
    )
    default_output = str(
        Path(settings["work_root"]) / settings["sample_id"] / "logs" / f"{stage}_%x_%j.out"
    )
    default_error = str(
        Path(settings["work_root"]) / settings["sample_id"] / "logs" / f"{stage}_%x_%j.err"
    )
    stage_settings["slurm_output"] = (
        pick(args.slurm_output, stage_slurm_cfg.get("output"))
        or stage_settings.get("slurm_output")
        or default_output
    )
    stage_settings["slurm_error"] = (
        pick(args.slurm_error, stage_slurm_cfg.get("error"))
        or stage_settings.get("slurm_error")
        or default_error
    )

    stage_settings["split_bams_slurm_partition"] = (
        pick(args.slurm_partition, split_bams_slurm_cfg.get("partition"))
        or stage_settings["slurm_partition"]
    )
    stage_settings["split_bams_slurm_mem"] = (
        pick(args.slurm_mem, split_bams_slurm_cfg.get("mem"))
        or stage_settings["slurm_mem"]
    )
    stage_settings["split_bams_slurm_cpus_per_task"] = (
        pick(args.slurm_cpus_per_task, split_bams_slurm_cfg.get("cpus_per_task"))
        or stage_settings["slurm_cpus_per_task"]
    )
    stage_settings["split_bams_slurm_output"] = (
        pick(args.slurm_output, split_bams_slurm_cfg.get("output"))
        or stage_settings["slurm_output"]
    )
    stage_settings["split_bams_slurm_error"] = (
        pick(args.slurm_error, split_bams_slurm_cfg.get("error"))
        or stage_settings["slurm_error"]
    )
    stage_settings["split_sort_slurm_partition"] = (
        pick(args.slurm_partition, split_sort_slurm_cfg.get("partition"))
        or stage_settings["slurm_partition"]
    )
    stage_settings["split_sort_slurm_mem"] = (
        pick(args.slurm_mem, split_sort_slurm_cfg.get("mem"))
        or stage_settings["slurm_mem"]
    )
    stage_settings["split_sort_slurm_cpus_per_task"] = (
        pick(args.slurm_cpus_per_task, split_sort_slurm_cfg.get("cpus_per_task"))
        or stage_settings["slurm_cpus_per_task"]
    )
    stage_settings["split_sort_slurm_output"] = (
        pick(args.slurm_output, split_sort_slurm_cfg.get("output"))
        or stage_settings["slurm_output"]
    )
    stage_settings["split_sort_slurm_error"] = (
        pick(args.slurm_error, split_sort_slurm_cfg.get("error"))
        or stage_settings["slurm_error"]
    )

    stage_settings["call_host_slurm_partition"] = (
        pick(args.slurm_partition, call_host_slurm_cfg.get("partition"))
        or stage_settings["slurm_partition"]
    )
    stage_settings["call_host_slurm_mem"] = (
        pick(args.slurm_mem, call_host_slurm_cfg.get("mem"))
        or stage_settings["slurm_mem"]
    )
    stage_settings["call_host_slurm_cpus_per_task"] = (
        pick(args.slurm_cpus_per_task, call_host_slurm_cfg.get("cpus_per_task"))
        or stage_settings["slurm_cpus_per_task"]
    )
    stage_settings["call_host_slurm_output"] = (
        pick(args.slurm_output, call_host_slurm_cfg.get("output"))
        or stage_settings["slurm_output"]
    )
    stage_settings["call_host_slurm_error"] = (
        pick(args.slurm_error, call_host_slurm_cfg.get("error"))
        or stage_settings["slurm_error"]
    )
    stage_settings["call_spike_slurm_partition"] = (
        pick(args.slurm_partition, call_spike_slurm_cfg.get("partition"))
        or stage_settings["slurm_partition"]
    )
    stage_settings["call_spike_slurm_mem"] = (
        pick(args.slurm_mem, call_spike_slurm_cfg.get("mem"))
        or stage_settings["slurm_mem"]
    )
    stage_settings["call_spike_slurm_cpus_per_task"] = (
        pick(args.slurm_cpus_per_task, call_spike_slurm_cfg.get("cpus_per_task"))
        or stage_settings["slurm_cpus_per_task"]
    )
    stage_settings["call_spike_slurm_output"] = (
        pick(args.slurm_output, call_spike_slurm_cfg.get("output"))
        or stage_settings["slurm_output"]
    )
    stage_settings["call_spike_slurm_error"] = (
        pick(args.slurm_error, call_spike_slurm_cfg.get("error"))
        or stage_settings["slurm_error"]
    )
    stage_settings["mbias_host_slurm_partition"] = (
        pick(args.slurm_partition, mbias_host_slurm_cfg.get("partition"))
        or stage_settings["slurm_partition"]
    )
    stage_settings["mbias_host_slurm_mem"] = (
        pick(args.slurm_mem, mbias_host_slurm_cfg.get("mem"))
        or stage_settings["slurm_mem"]
    )
    stage_settings["mbias_host_slurm_cpus_per_task"] = (
        pick(args.slurm_cpus_per_task, mbias_host_slurm_cfg.get("cpus_per_task"))
        or stage_settings["slurm_cpus_per_task"]
    )
    stage_settings["mbias_host_slurm_output"] = (
        pick(args.slurm_output, mbias_host_slurm_cfg.get("output"))
        or stage_settings["slurm_output"]
    )
    stage_settings["mbias_host_slurm_error"] = (
        pick(args.slurm_error, mbias_host_slurm_cfg.get("error"))
        or stage_settings["slurm_error"]
    )
    stage_settings["mbias_spike_slurm_partition"] = (
        pick(args.slurm_partition, mbias_spike_slurm_cfg.get("partition"))
        or stage_settings["slurm_partition"]
    )
    stage_settings["mbias_spike_slurm_mem"] = (
        pick(args.slurm_mem, mbias_spike_slurm_cfg.get("mem"))
        or stage_settings["slurm_mem"]
    )
    stage_settings["mbias_spike_slurm_cpus_per_task"] = (
        pick(args.slurm_cpus_per_task, mbias_spike_slurm_cfg.get("cpus_per_task"))
        or stage_settings["slurm_cpus_per_task"]
    )
    stage_settings["mbias_spike_slurm_output"] = (
        pick(args.slurm_output, mbias_spike_slurm_cfg.get("output"))
        or stage_settings["slurm_output"]
    )
    stage_settings["mbias_spike_slurm_error"] = (
        pick(args.slurm_error, mbias_spike_slurm_cfg.get("error"))
        or stage_settings["slurm_error"]
    )
    return stage_settings


def resolve_step_slurm_cfg(
    stage_slurm_cfg: dict,
    step_name: str,
    step_keys: set[str],
) -> dict:
    if not isinstance(stage_slurm_cfg, dict):
        raise ValueError("selected slurm config must be an object")
    if not any(key in stage_slurm_cfg for key in step_keys):
        return stage_slurm_cfg
    step_cfg = stage_slurm_cfg.get(step_name, {})
    if step_cfg is None:
        step_cfg = {}
    if not isinstance(step_cfg, dict):
        raise ValueError(f"slurm split config '{step_name}' must be an object")
    base_cfg = {
        key: value for key, value in stage_slurm_cfg.items() if key not in step_keys
    }
    return {**base_cfg, **step_cfg}


def resolve_settings(args: argparse.Namespace) -> dict:
    cfg: dict = {}
    if args.workflow_config:
        cfg = load_workflow_config(args.workflow_config)

    slurm_cfg_raw = cfg.get("slurm", {})
    if slurm_cfg_raw is None:
        slurm_cfg_raw = {}
    if not isinstance(slurm_cfg_raw, dict):
        raise ValueError("workflow config key 'slurm' must be an object")
    stage = pick(args.stage, cfg.get("stage")) or "fastp_split"
    if stage == "all":
        slurm_stage_key = STAGE_SEQUENCE[0]
    elif stage in STAGE_SEQUENCE or stage in ("aggregate", *METHSCAN_STAGES):
        slurm_stage_key = stage
    else:
        slurm_stage_key = STAGE_SEQUENCE[0]
    stage_slurm_cfg = select_stage_slurm_cfg(slurm_cfg_raw, slurm_stage_key)
    split_bams_slurm_cfg = resolve_step_slurm_cfg(
        stage_slurm_cfg,
        "split_bams",
        {"split_bams", "sort"},
    )
    split_sort_slurm_cfg = resolve_step_slurm_cfg(
        stage_slurm_cfg,
        "sort",
        {"split_bams", "sort"},
    )
    call_host_slurm_cfg = resolve_step_slurm_cfg(
        stage_slurm_cfg,
        "host",
        {"host", "spike"},
    )
    call_spike_slurm_cfg = resolve_step_slurm_cfg(
        stage_slurm_cfg,
        "spike",
        {"host", "spike"},
    )
    mbias_host_slurm_cfg = resolve_step_slurm_cfg(
        stage_slurm_cfg,
        "host",
        {"host", "spike"},
    )
    mbias_spike_slurm_cfg = resolve_step_slurm_cfg(
        stage_slurm_cfg,
        "spike",
        {"host", "spike"},
    )

    settings = {
        "runner": pick(args.runner, cfg.get("runner")),
        "stage": stage,
        "sample_id": pick(args.sample_id, cfg.get("sample_id")),
        "r1": pick(args.r1, cfg.get("r1")),
        "r2": pick(args.r2, cfg.get("r2")),
        "work_root": pick(args.work_root, cfg.get("work_root")),
        "fastp_threads": pick(args.fastp_threads, cfg.get("fastp_threads")),
        "number_of_split_parts": pick(
            args.number_of_split_parts, cfg.get("number_of_split_parts")
        ),
        "fastp_bin": pick(args.fastp_bin, cfg.get("fastp_bin")),
        "barcode1_whitelist": pick(
            args.barcode1_whitelist, cfg.get("barcode1_whitelist")
        ),
        "barcode2_whitelist": pick(
            args.barcode2_whitelist, cfg.get("barcode2_whitelist")
        ),
        "linker_bc": pick(
            args.linker_bc,
            cfg.get("linker_bc", cfg.get("linker2")),
        ),
        "insert_left": pick(
            args.insert_left,
            cfg.get("insert_left", cfg.get("tn5")),
        ),
        "linker1": pick(args.linker1, cfg.get("linker1")),
        "linker2": pick(args.linker2, cfg.get("linker2")),
        "tn5": pick(args.tn5, cfg.get("tn5")),
        "linker_edit_distance": pick(
            args.linker_edit_distance, cfg.get("linker_edit_distance")
        ),
        "barcode_hamming_distance": pick(
            args.barcode_hamming_distance, cfg.get("barcode_hamming_distance")
        ),
        "gzip_level": pick(args.gzip_level, cfg.get("gzip_level")),
        "spike_edit_distance": pick(
            args.spike_edit_distance, cfg.get("spike_edit_distance")
        ),
        "require_c_all_t": pick(args.require_c_all_t, cfg.get("require_c_all_t")),
        "bwa_index": pick(args.bwa_index, cfg.get("bwa_index")),
        "bwa_threads": pick(args.bwa_threads, cfg.get("bwa_threads")),
        "bwa_bin": pick(args.bwa_bin, cfg.get("bwa_bin")),
        "sinto_bin": pick(args.sinto_bin, cfg.get("sinto_bin")),
        "samtools_bin": pick(args.samtools_bin, cfg.get("samtools_bin")),
        "samtools_threads": pick(args.samtools_threads, cfg.get("samtools_threads")),
        "host_sort_mem": pick(args.host_sort_mem, cfg.get("host_sort_mem")),
        "split_barcodes": pick(
            args.split_barcodes,
            cfg.get("split_barcodes", cfg.get("barcode1_whitelist")),
        ),
        "split_cb_tag": pick(args.split_cb_tag, cfg.get("split_cb_tag")),
        "split_threads_read": pick(
            args.split_threads_read, cfg.get("split_threads_read")
        ),
        "split_threads_write": pick(
            args.split_threads_write, cfg.get("split_threads_write")
        ),
        "split_smoke": pick(args.split_smoke, cfg.get("split_smoke")),
        "split_sort_jobs": pick(args.split_sort_jobs, cfg.get("split_sort_jobs")),
        "mbias_host_subsample_fraction": pick(
            args.mbias_host_subsample_fraction, cfg.get("mbias_host_subsample_fraction")
        ),
        "mbias_mode": pick(args.mbias_mode, cfg.get("mbias_mode")),
        "mbias_max_cycle": pick(args.mbias_max_cycle, cfg.get("mbias_max_cycle")),
        "mbias_min_mapping_quality": pick(
            args.mbias_min_mapping_quality, cfg.get("mbias_min_mapping_quality")
        ),
        "mbias_script": pick(args.mbias_script, cfg.get("mbias_script")),
        "call_reference_file": pick(
            args.call_reference_file,
            cfg.get("call_reference_file", cfg.get("bwa_index")),
        ),
        "call_chromosomes": pick(args.call_chromosomes, cfg.get("call_chromosomes")),
        "call_mito_chromosomes": pick(
            args.call_mito_chromosomes, cfg.get("call_mito_chromosomes")
        ),
        "call_jobs": pick(args.call_jobs, cfg.get("call_jobs")),
        "call_min_base_quality": pick(
            args.call_min_base_quality, cfg.get("call_min_base_quality")
        ),
        "call_min_mapping_quality": pick(
            args.call_min_mapping_quality, cfg.get("call_min_mapping_quality")
        ),
        "call_sample_size": pick(args.call_sample_size, cfg.get("call_sample_size")),
        "call_max_depth": pick(args.call_max_depth, cfg.get("call_max_depth")),
        "call_batch_size": pick(args.call_batch_size, cfg.get("call_batch_size")),
        "call_mode": pick(args.call_mode, cfg.get("call_mode")),
        "call_context_mode": pick(
            args.call_context_mode,
            cfg.get("call_context_mode"),
        ),
        "call_r1_left_trimming": pick(
            args.call_r1_left_trimming, cfg.get("call_r1_left_trimming")
        ),
        "call_r1_right_trimming": pick(
            args.call_r1_right_trimming, cfg.get("call_r1_right_trimming")
        ),
        "call_r2_left_trimming": pick(
            args.call_r2_left_trimming, cfg.get("call_r2_left_trimming")
        ),
        "call_r2_right_trimming": pick(
            args.call_r2_right_trimming, cfg.get("call_r2_right_trimming")
        ),
        "call_caller_script": pick(args.call_caller_script, cfg.get("call_caller_script")),
        "call_ch_caller_script": pick(
            args.call_ch_caller_script,
            cfg.get("call_ch_caller_script"),
        ),
        "saturation_script": pick(args.saturation_script, cfg.get("saturation_script")),
        "saturation_reads_threshold": pick(
            args.saturation_reads_threshold,
            cfg.get("saturation_reads_threshold"),
        ),
        "summary_script": pick(args.summary_script, cfg.get("summary_script")),
        "aggregate_script": pick(args.aggregate_script, cfg.get("aggregate_script")),
        "aggregate_sort_mem": pick(
            args.aggregate_sort_mem,
            cfg.get("aggregate_sort_mem", "8G"),
        ),
        "methscan_run_script": pick(args.methscan_run_script, cfg.get("methscan_run_script")),
        "methscan_prepare_chunksize": pick(
            args.methscan_prepare_chunksize,
            cfg.get("methscan_prepare_chunksize"),
        ),
        "methscan_filter_min_sites": pick(
            args.methscan_filter_min_sites,
            cfg.get("methscan_filter_min_sites"),
        ),
        "methscan_tss_bed": pick(args.methscan_tss_bed, cfg.get("methscan_tss_bed")),
        "methscan_profile_strand_column": pick(
            args.methscan_profile_strand_column,
            cfg.get("methscan_profile_strand_column"),
        ),
        "methscan_profile_prepared_subdir": pick(
            args.methscan_profile_prepared_subdir,
            cfg.get("methscan_profile_prepared_subdir"),
        ),
        "methscan_profile_csv": pick(args.methscan_profile_csv, cfg.get("methscan_profile_csv")),
        "methscan_smooth_bandwidth": pick(
            args.methscan_smooth_bandwidth,
            cfg.get("methscan_smooth_bandwidth"),
        ),
        "methscan_smooth_use_weights": pick(
            args.methscan_smooth_use_weights,
            cfg.get("methscan_smooth_use_weights"),
        ),
        "methscan_scan_threads": pick(
            args.methscan_scan_threads,
            cfg.get("methscan_scan_threads"),
        ),
        "methscan_vmrs_bed": pick(args.methscan_vmrs_bed, cfg.get("methscan_vmrs_bed")),
        "methscan_matrix_threads": pick(
            args.methscan_matrix_threads,
            cfg.get("methscan_matrix_threads"),
        ),
        "methscan_matrix_sparse": pick(
            args.methscan_matrix_sparse,
            cfg.get("methscan_matrix_sparse"),
        ),
        "methscan_matrix_prefix": pick(
            args.methscan_matrix_prefix,
            cfg.get("methscan_matrix_prefix"),
        ),
        "methscan_pixi_manifest": pick(
            args.methscan_pixi_manifest,
            cfg.get("methscan_pixi_manifest"),
        ),
        "spike_in_index": normalize_spike_in_index(
            pick(args.spike_in_index, cfg.get("spike_in_index"))
        ),
        "slurm_cfg_raw": slurm_cfg_raw,
        "slurm_partition": pick(args.slurm_partition, stage_slurm_cfg.get("partition")),
        "slurm_mem": pick(args.slurm_mem, stage_slurm_cfg.get("mem")),
        "slurm_cpus_per_task": pick(
            args.slurm_cpus_per_task, stage_slurm_cfg.get("cpus_per_task")
        ),
        "slurm_output": pick(args.slurm_output, stage_slurm_cfg.get("output")),
        "slurm_error": pick(args.slurm_error, stage_slurm_cfg.get("error")),
        "split_bams_slurm_partition": pick(
            args.slurm_partition, split_bams_slurm_cfg.get("partition")
        ),
        "split_bams_slurm_mem": pick(args.slurm_mem, split_bams_slurm_cfg.get("mem")),
        "split_bams_slurm_cpus_per_task": pick(
            args.slurm_cpus_per_task, split_bams_slurm_cfg.get("cpus_per_task")
        ),
        "split_bams_slurm_output": pick(
            args.slurm_output, split_bams_slurm_cfg.get("output")
        ),
        "split_bams_slurm_error": pick(
            args.slurm_error, split_bams_slurm_cfg.get("error")
        ),
        "split_sort_slurm_partition": pick(
            args.slurm_partition, split_sort_slurm_cfg.get("partition")
        ),
        "split_sort_slurm_mem": pick(args.slurm_mem, split_sort_slurm_cfg.get("mem")),
        "split_sort_slurm_cpus_per_task": pick(
            args.slurm_cpus_per_task, split_sort_slurm_cfg.get("cpus_per_task")
        ),
        "split_sort_slurm_output": pick(
            args.slurm_output, split_sort_slurm_cfg.get("output")
        ),
        "split_sort_slurm_error": pick(
            args.slurm_error, split_sort_slurm_cfg.get("error")
        ),
        "call_host_slurm_partition": pick(
            args.slurm_partition, call_host_slurm_cfg.get("partition")
        ),
        "call_host_slurm_mem": pick(args.slurm_mem, call_host_slurm_cfg.get("mem")),
        "call_host_slurm_cpus_per_task": pick(
            args.slurm_cpus_per_task, call_host_slurm_cfg.get("cpus_per_task")
        ),
        "call_host_slurm_output": pick(
            args.slurm_output, call_host_slurm_cfg.get("output")
        ),
        "call_host_slurm_error": pick(
            args.slurm_error, call_host_slurm_cfg.get("error")
        ),
        "call_spike_slurm_partition": pick(
            args.slurm_partition, call_spike_slurm_cfg.get("partition")
        ),
        "call_spike_slurm_mem": pick(args.slurm_mem, call_spike_slurm_cfg.get("mem")),
        "call_spike_slurm_cpus_per_task": pick(
            args.slurm_cpus_per_task, call_spike_slurm_cfg.get("cpus_per_task")
        ),
        "call_spike_slurm_output": pick(
            args.slurm_output, call_spike_slurm_cfg.get("output")
        ),
        "call_spike_slurm_error": pick(
            args.slurm_error, call_spike_slurm_cfg.get("error")
        ),
        "mbias_host_slurm_partition": pick(
            args.slurm_partition, mbias_host_slurm_cfg.get("partition")
        ),
        "mbias_host_slurm_mem": pick(args.slurm_mem, mbias_host_slurm_cfg.get("mem")),
        "mbias_host_slurm_cpus_per_task": pick(
            args.slurm_cpus_per_task, mbias_host_slurm_cfg.get("cpus_per_task")
        ),
        "mbias_host_slurm_output": pick(
            args.slurm_output, mbias_host_slurm_cfg.get("output")
        ),
        "mbias_host_slurm_error": pick(
            args.slurm_error, mbias_host_slurm_cfg.get("error")
        ),
        "mbias_spike_slurm_partition": pick(
            args.slurm_partition, mbias_spike_slurm_cfg.get("partition")
        ),
        "mbias_spike_slurm_mem": pick(args.slurm_mem, mbias_spike_slurm_cfg.get("mem")),
        "mbias_spike_slurm_cpus_per_task": pick(
            args.slurm_cpus_per_task, mbias_spike_slurm_cfg.get("cpus_per_task")
        ),
        "mbias_spike_slurm_output": pick(
            args.slurm_output, mbias_spike_slurm_cfg.get("output")
        ),
        "mbias_spike_slurm_error": pick(
            args.slurm_error, mbias_spike_slurm_cfg.get("error")
        ),
        "submit": args.submit,
        "dry_run": args.dry_run,
    }

    if stage not in STAGE_CHOICES:
        raise ValueError(f"unsupported stage: {stage}")

    settings["work_root"] = settings["work_root"] or "work"
    settings["fastp_threads"] = settings["fastp_threads"] or 8
    if settings["number_of_split_parts"] is not None:
        settings["number_of_split_parts"] = int(settings["number_of_split_parts"])
        if settings["number_of_split_parts"] <= 0:
            raise ValueError("number_of_split_parts must be > 0")
    settings["fastp_bin"] = normalize_executable_setting(
        settings["fastp_bin"], "fastp"
    )
    settings["linker1"] = settings["linker1"] or "GTGGCCGATGTTTCG"
    settings["linker2"] = (
        settings["linker2"] or "ATCCACGTGCTTGAGAGGCCAGAGCATTCG"
    )
    settings["tn5"] = settings["tn5"] or "CATCGGCGTACGACTAGATGTGTATAAGAGACAG"
    settings["linker_bc"] = (
        settings["linker_bc"]
        or settings["linker2"]
        or "ATCCACGTGCTTGAGAGGCCAGAGCATTCG"
    )
    settings["insert_left"] = (
        settings["insert_left"]
        or settings["tn5"]
        or "CATCGGCGTACGACTAGATGTGTATAAGAGACAG"
    )
    settings["linker_edit_distance"] = (
        int(settings["linker_edit_distance"])
        if settings["linker_edit_distance"] is not None
        else 1
    )
    settings["barcode_hamming_distance"] = (
        int(settings["barcode_hamming_distance"])
        if settings["barcode_hamming_distance"] is not None
        else 1
    )
    if settings["linker_edit_distance"] < 0:
        raise ValueError("linker_edit_distance must be >= 0")
    if settings["barcode_hamming_distance"] < 0:
        raise ValueError("barcode_hamming_distance must be >= 0")
    settings["gzip_level"] = (
        int(settings["gzip_level"]) if settings["gzip_level"] is not None else 6
    )
    if settings["gzip_level"] < 0 or settings["gzip_level"] > 9:
        raise ValueError("gzip_level must be between 0 and 9")
    settings["spike_edit_distance"] = (
        int(settings["spike_edit_distance"])
        if settings["spike_edit_distance"] is not None
        else 1
    )
    if settings["spike_edit_distance"] < 0:
        raise ValueError("spike_edit_distance must be >= 0")
    settings["require_c_all_t"] = (
        "" if settings["require_c_all_t"] is None else str(settings["require_c_all_t"])
    )
    settings["bwa_threads"] = (
        int(settings["bwa_threads"]) if settings["bwa_threads"] is not None else 2
    )
    if settings["bwa_threads"] <= 0:
        raise ValueError("bwa_threads must be > 0")
    settings["bwa_bin"] = normalize_executable_setting(settings["bwa_bin"], "bwa")
    settings["sinto_bin"] = normalize_executable_setting(
        settings["sinto_bin"], "sinto"
    )
    settings["samtools_bin"] = normalize_executable_setting(
        settings["samtools_bin"], "samtools"
    )
    settings["samtools_threads"] = (
        int(settings["samtools_threads"])
        if settings["samtools_threads"] is not None
        else 4
    )
    if settings["samtools_threads"] <= 0:
        raise ValueError("samtools_threads must be > 0")
    settings["host_sort_mem"] = settings["host_sort_mem"] or "16G"
    settings["split_barcodes"] = (
        settings["split_barcodes"] or settings["barcode1_whitelist"]
    )
    settings["split_cb_tag"] = settings["split_cb_tag"] or "CB"
    settings["split_threads_read"] = (
        int(settings["split_threads_read"])
        if settings["split_threads_read"] is not None
        else 1
    )
    settings["split_threads_write"] = (
        int(settings["split_threads_write"])
        if settings["split_threads_write"] is not None
        else 1
    )
    if settings["split_threads_read"] <= 0:
        raise ValueError("split_threads_read must be > 0")
    if settings["split_threads_write"] <= 0:
        raise ValueError("split_threads_write must be > 0")
    settings["split_smoke"] = bool(settings["split_smoke"])
    settings["split_sort_jobs"] = (
        int(settings["split_sort_jobs"]) if settings["split_sort_jobs"] is not None else 8
    )
    if settings["split_sort_jobs"] <= 0:
        raise ValueError("split_sort_jobs must be > 0")
    settings["mbias_host_subsample_fraction"] = (
        float(settings["mbias_host_subsample_fraction"])
        if settings["mbias_host_subsample_fraction"] is not None
        else 0.1
    )
    if (
        settings["mbias_host_subsample_fraction"] <= 0
        or settings["mbias_host_subsample_fraction"] > 1
    ):
        raise ValueError("mbias_host_subsample_fraction must be in (0, 1]")
    settings["mbias_mode"] = settings["mbias_mode"] or "spike"
    if settings["mbias_mode"] not in {"all", "host", "spike"}:
        raise ValueError("mbias_mode must be one of: all, host, spike")
    if settings["mbias_mode"] in {"all", "host"} and not settings.get("call_chromosomes"):
        raise ValueError(
            "mbias stage requires call_chromosomes when mbias_mode is all or host"
        )
    settings["mbias_max_cycle"] = (
        int(settings["mbias_max_cycle"]) if settings["mbias_max_cycle"] is not None else 150
    )
    if settings["mbias_max_cycle"] <= 0:
        raise ValueError("mbias_max_cycle must be > 0")
    settings["mbias_min_mapping_quality"] = (
        int(settings["mbias_min_mapping_quality"])
        if settings["mbias_min_mapping_quality"] is not None
        else 1
    )
    if settings["mbias_min_mapping_quality"] < 0:
        raise ValueError("mbias_min_mapping_quality must be >= 0")
    settings["mbias_script"] = settings["mbias_script"] or "scripts/mbias.py"
    settings["call_mito_chromosomes"] = settings["call_mito_chromosomes"] or "chrM"
    settings["call_jobs"] = (
        int(settings["call_jobs"]) if settings["call_jobs"] is not None else 8
    )
    if settings["call_jobs"] <= 0:
        raise ValueError("call_jobs must be > 0")
    settings["call_min_base_quality"] = (
        int(settings["call_min_base_quality"])
        if settings["call_min_base_quality"] is not None
        else 30
    )
    if settings["call_min_base_quality"] < 0:
        raise ValueError("call_min_base_quality must be >= 0")
    settings["call_min_mapping_quality"] = (
        int(settings["call_min_mapping_quality"])
        if settings["call_min_mapping_quality"] is not None
        else 1
    )
    if settings["call_min_mapping_quality"] < 0:
        raise ValueError("call_min_mapping_quality must be >= 0")
    if settings["call_sample_size"] is not None:
        settings["call_sample_size"] = int(settings["call_sample_size"])
        if settings["call_sample_size"] <= 0:
            raise ValueError("call_sample_size must be > 0 when provided")
    settings["call_max_depth"] = (
        int(settings["call_max_depth"]) if settings["call_max_depth"] is not None else 250
    )
    if settings["call_max_depth"] <= 0:
        raise ValueError("call_max_depth must be > 0")
    settings["call_batch_size"] = (
        int(settings["call_batch_size"])
        if settings["call_batch_size"] is not None
        else 10_000_000
    )
    if settings["call_batch_size"] <= 0:
        raise ValueError("call_batch_size must be > 0")
    settings["call_mode"] = settings["call_mode"] or "all"
    if settings["call_mode"] not in {"all", "host", "spike"}:
        raise ValueError("call_mode must be one of: all, host, spike")
    settings["call_context_mode"] = settings["call_context_mode"] or "cg"
    if settings["call_context_mode"] not in {"cg", "ch", "both"}:
        raise ValueError("call_context_mode must be one of: cg, ch, both")
    settings["call_r1_left_trimming"] = (
        int(settings["call_r1_left_trimming"])
        if settings["call_r1_left_trimming"] is not None
        else 0
    )
    if settings["call_r1_left_trimming"] < 0:
        raise ValueError("call_r1_left_trimming must be >= 0")
    settings["call_r1_right_trimming"] = (
        int(settings["call_r1_right_trimming"])
        if settings["call_r1_right_trimming"] is not None
        else 0
    )
    if settings["call_r1_right_trimming"] < 0:
        raise ValueError("call_r1_right_trimming must be >= 0")
    settings["call_r2_left_trimming"] = (
        int(settings["call_r2_left_trimming"])
        if settings["call_r2_left_trimming"] is not None
        else 0
    )
    if settings["call_r2_left_trimming"] < 0:
        raise ValueError("call_r2_left_trimming must be >= 0")
    settings["call_r2_right_trimming"] = (
        int(settings["call_r2_right_trimming"])
        if settings["call_r2_right_trimming"] is not None
        else 0
    )
    if settings["call_r2_right_trimming"] < 0:
        raise ValueError("call_r2_right_trimming must be >= 0")
    call_trimming_cfg = cfg.get("call_trimming")
    if call_trimming_cfg is None:
        call_trimming_cfg = {}
    if not isinstance(call_trimming_cfg, dict):
        raise ValueError("call_trimming must be an object")
    flat_call_trim = {
        "r1_left": settings["call_r1_left_trimming"],
        "r1_right": settings["call_r1_right_trimming"],
        "r2_left": settings["call_r2_left_trimming"],
        "r2_right": settings["call_r2_right_trimming"],
    }
    settings["call_trim_host"] = resolve_call_trim_target(
        call_trimming_cfg.get("host"), flat_call_trim, "host"
    )
    settings["call_trim_spike"] = resolve_call_trim_target(
        call_trimming_cfg.get("spike"), flat_call_trim, "spike"
    )
    settings["call_caller_script"] = (
        settings["call_caller_script"] or "scripts/methy_caller.py"
    )
    settings["call_ch_caller_script"] = (
        settings["call_ch_caller_script"] or "scripts/methy_caller_CH.py"
    )
    settings["saturation_script"] = settings["saturation_script"] or "scripts/saturation.py"
    settings["saturation_reads_threshold"] = (
        float(settings["saturation_reads_threshold"])
        if settings["saturation_reads_threshold"] is not None
        else 1_000_000.0
    )
    if settings["saturation_reads_threshold"] <= 0:
        raise ValueError("saturation_reads_threshold must be > 0")
    settings["summary_script"] = settings["summary_script"] or "scripts-v2/summary.py"
    settings["aggregate_script"] = (
        settings["aggregate_script"] or "scripts/aggregate.py"
    )
    settings["methscan_run_script"] = settings.get("methscan_run_script") or "scripts/methscan_run.py"
    settings["methscan_prepare_chunksize"] = (
        int(settings["methscan_prepare_chunksize"])
        if settings.get("methscan_prepare_chunksize") is not None
        else 10_000_000
    )
    if settings["methscan_prepare_chunksize"] < 1:
        raise ValueError("methscan_prepare_chunksize must be >= 1")
    settings["methscan_filter_min_sites"] = (
        int(settings["methscan_filter_min_sites"])
        if settings.get("methscan_filter_min_sites") is not None
        else 50_000
    )
    if settings["methscan_filter_min_sites"] < 1:
        raise ValueError("methscan_filter_min_sites must be >= 1")
    settings["methscan_profile_strand_column"] = (
        int(settings["methscan_profile_strand_column"])
        if settings.get("methscan_profile_strand_column") is not None
        else 6
    )
    if settings["methscan_profile_strand_column"] < 1:
        raise ValueError("methscan_profile_strand_column must be >= 1")
    settings["methscan_profile_prepared_subdir"] = (
        settings.get("methscan_profile_prepared_subdir") or "compact"
    )
    if settings["methscan_profile_prepared_subdir"] not in ("compact", "filter"):
        raise ValueError(
            "methscan_profile_prepared_subdir must be 'compact' or 'filter'"
        )
    settings["methscan_scan_threads"] = (
        int(settings["methscan_scan_threads"])
        if settings.get("methscan_scan_threads") is not None
        else 10
    )
    if settings["methscan_scan_threads"] < 1:
        raise ValueError("methscan_scan_threads must be >= 1")
    settings["methscan_matrix_threads"] = (
        int(settings["methscan_matrix_threads"])
        if settings.get("methscan_matrix_threads") is not None
        else 10
    )
    if settings["methscan_matrix_threads"] < 1:
        raise ValueError("methscan_matrix_threads must be >= 1")
    settings["methscan_smooth_use_weights"] = bool(
        settings.get("methscan_smooth_use_weights")
    )
    settings["methscan_matrix_sparse"] = bool(settings.get("methscan_matrix_sparse"))
    settings["slurm_partition"] = settings["slurm_partition"] or "cpu"
    settings["slurm_mem"] = settings["slurm_mem"] or "16G"
    settings["slurm_cpus_per_task"] = settings["slurm_cpus_per_task"] or 8
    settings["slurm_output"] = settings["slurm_output"] or str(
        Path(settings["work_root"])
        / settings["sample_id"]
        / "logs"
        / f"{stage}_%x_%j.out"
    )
    settings["slurm_error"] = settings["slurm_error"] or str(
        Path(settings["work_root"])
        / settings["sample_id"]
        / "logs"
        / f"{stage}_%x_%j.err"
    )
    settings["split_bams_slurm_partition"] = (
        settings["split_bams_slurm_partition"] or settings["slurm_partition"]
    )
    settings["split_bams_slurm_mem"] = (
        settings["split_bams_slurm_mem"] or settings["slurm_mem"]
    )
    settings["split_bams_slurm_cpus_per_task"] = (
        settings["split_bams_slurm_cpus_per_task"] or settings["slurm_cpus_per_task"]
    )
    settings["split_bams_slurm_output"] = (
        settings["split_bams_slurm_output"] or settings["slurm_output"]
    )
    settings["split_bams_slurm_error"] = (
        settings["split_bams_slurm_error"] or settings["slurm_error"]
    )
    settings["split_sort_slurm_partition"] = (
        settings["split_sort_slurm_partition"] or settings["slurm_partition"]
    )
    settings["split_sort_slurm_mem"] = (
        settings["split_sort_slurm_mem"] or settings["slurm_mem"]
    )
    settings["split_sort_slurm_cpus_per_task"] = (
        settings["split_sort_slurm_cpus_per_task"] or settings["slurm_cpus_per_task"]
    )
    settings["split_sort_slurm_output"] = (
        settings["split_sort_slurm_output"] or settings["slurm_output"]
    )
    settings["split_sort_slurm_error"] = (
        settings["split_sort_slurm_error"] or settings["slurm_error"]
    )
    settings["call_host_slurm_partition"] = (
        settings["call_host_slurm_partition"] or settings["slurm_partition"]
    )
    settings["call_host_slurm_mem"] = (
        settings["call_host_slurm_mem"] or settings["slurm_mem"]
    )
    settings["call_host_slurm_cpus_per_task"] = (
        settings["call_host_slurm_cpus_per_task"] or settings["slurm_cpus_per_task"]
    )
    settings["call_host_slurm_output"] = (
        settings["call_host_slurm_output"] or settings["slurm_output"]
    )
    settings["call_host_slurm_error"] = (
        settings["call_host_slurm_error"] or settings["slurm_error"]
    )
    settings["call_spike_slurm_partition"] = (
        settings["call_spike_slurm_partition"] or settings["slurm_partition"]
    )
    settings["call_spike_slurm_mem"] = (
        settings["call_spike_slurm_mem"] or settings["slurm_mem"]
    )
    settings["call_spike_slurm_cpus_per_task"] = (
        settings["call_spike_slurm_cpus_per_task"] or settings["slurm_cpus_per_task"]
    )
    settings["call_spike_slurm_output"] = (
        settings["call_spike_slurm_output"] or settings["slurm_output"]
    )
    settings["call_spike_slurm_error"] = (
        settings["call_spike_slurm_error"] or settings["slurm_error"]
    )
    settings["mbias_host_slurm_partition"] = (
        settings["mbias_host_slurm_partition"] or settings["slurm_partition"]
    )
    settings["mbias_host_slurm_mem"] = (
        settings["mbias_host_slurm_mem"] or settings["slurm_mem"]
    )
    settings["mbias_host_slurm_cpus_per_task"] = (
        settings["mbias_host_slurm_cpus_per_task"] or settings["slurm_cpus_per_task"]
    )
    settings["mbias_host_slurm_output"] = (
        settings["mbias_host_slurm_output"] or settings["slurm_output"]
    )
    settings["mbias_host_slurm_error"] = (
        settings["mbias_host_slurm_error"] or settings["slurm_error"]
    )
    settings["mbias_spike_slurm_partition"] = (
        settings["mbias_spike_slurm_partition"] or settings["slurm_partition"]
    )
    settings["mbias_spike_slurm_mem"] = (
        settings["mbias_spike_slurm_mem"] or settings["slurm_mem"]
    )
    settings["mbias_spike_slurm_cpus_per_task"] = (
        settings["mbias_spike_slurm_cpus_per_task"] or settings["slurm_cpus_per_task"]
    )
    settings["mbias_spike_slurm_output"] = (
        settings["mbias_spike_slurm_output"] or settings["slurm_output"]
    )
    settings["mbias_spike_slurm_error"] = (
        settings["mbias_spike_slurm_error"] or settings["slurm_error"]
    )

    if stage == "all":
        for stage_name in STAGE_SEQUENCE:
            validate_required_for_stage(stage_name, settings)
    else:
        validate_required_for_stage(stage, settings)

    return settings


def generate_local_script(command: str, output_path: Path) -> None:
    content = (
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n\n"
        f"{command}\n"
    )
    write_text(output_path, content)
    output_path.chmod(0o755)


def generate_slurm_script(
    command: str, output_path: Path, log_dir: Path, args: argparse.Namespace
) -> None:
    log_dir.mkdir(parents=True, exist_ok=True)
    lines = [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name={args.job_name}",
        f"#SBATCH --partition={args.slurm_partition}",
        f"#SBATCH --cpus-per-task={args.slurm_cpus_per_task}",
        f"#SBATCH --mem={args.slurm_mem}",
        f"#SBATCH --output={args.slurm_output}",
        f"#SBATCH --error={args.slurm_error}",
        "",
        "set -euo pipefail",
        "",
    ]
    if args.module_line:
        lines.extend([args.module_line, ""])
    lines.append(command)
    lines.append("")
    content = "\n".join(lines)
    write_text(output_path, content)
    output_path.chmod(0o755)


def submit_script(path: Path, runner: str) -> None:
    if runner == "local":
        subprocess.run(["bash", str(path)], check=True)
    else:
        subprocess.run(["sbatch", str(path)], check=True)


def submit_slurm_script(path: Path, dependency_job_id: str | None = None) -> str:
    command = ["sbatch"]
    if dependency_job_id:
        command.append(f"--dependency=afterok:{dependency_job_id}")
    command.append(str(path))
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    stdout = completed.stdout.strip()
    if stdout:
        print(stdout)
    tokens = stdout.split()
    if not tokens:
        raise ValueError(f"failed to parse sbatch output for: {path}")
    return tokens[-1]


def parse_generated_paths(command_output: str) -> list[Path]:
    generated: list[Path] = []
    for line in command_output.splitlines():
        prefix = "[make_cmd] generated="
        if line.startswith(prefix):
            generated.append(Path(line[len(prefix) :].strip()))
    return generated


def join_dependency_ids(job_ids: list[str]) -> str:
    return ":".join(job_id for job_id in job_ids if job_id)


def build_stage_passthrough_args(argv: list[str]) -> list[str]:
    passthrough: list[str] = []
    flags_without_value = {
        "--split-smoke",
        "--submit",
        "--dry-run",
        "--skip-workdir-input-checks",
    }
    index = 0
    while index < len(argv):
        token = argv[index]
        if token in {"--stage", "--runner"}:
            index += 2
            continue
        if token.startswith("--stage=") or token.startswith("--runner="):
            index += 1
            continue
        if token in {"--submit", "--dry-run", "--skip-workdir-input-checks"}:
            index += 1
            continue
        if token in flags_without_value:
            passthrough.append(token)
            index += 1
            continue
        if token.startswith("--"):
            passthrough.append(token)
            if index + 1 < len(argv):
                passthrough.append(argv[index + 1])
            index += 2
            continue
        passthrough.append(token)
        index += 1
    return passthrough


def generate_local_driver_script(stage_scripts: list[tuple[str, list[Path]]], output_path: Path) -> None:
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        'SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
        "",
    ]
    for _, scripts in stage_scripts:
        for script_path in scripts:
            lines.append(f'bash "$SCRIPT_DIR/{script_path.name}"')
    lines.append("")
    write_text(output_path, "\n".join(lines))
    output_path.chmod(0o755)


def generate_split_slurm_submit_script(
    split_script_path: Path, sort_script_path: Path, output_path: Path
) -> None:
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        'SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
        "",
        'split_out="$(sbatch "$SCRIPT_DIR/' + split_script_path.name + '")"',
        'echo "$split_out"',
        'jid_split_bams="${split_out##* }"',
        "",
        'sort_out="$(sbatch --dependency=afterok:${jid_split_bams} "$SCRIPT_DIR/'
        + sort_script_path.name
        + '")"',
        'echo "$sort_out"',
        'jid_split_sort="${sort_out##* }"',
        "",
        'echo "[split.submit] done split_bams=${jid_split_bams} split_sort=${jid_split_sort}"',
        "",
    ]
    write_text(output_path, "\n".join(lines))
    output_path.chmod(0o755)


def generate_slurm_driver_script(
    stage_scripts: list[tuple[str, list[Path]]],
    output_path: Path,
    log_dir: Path,
    settings: dict,
) -> None:
    lines = [
        "submit_with_dep() {",
        '  local script_path="$1"',
        '  local dep_chain="$2"',
        "  local out",
        '  if [[ -n "$dep_chain" ]]; then',
        '    out="$(sbatch --dependency=afterok:${dep_chain} "$script_path")"',
        "  else",
        '    out="$(sbatch "$script_path")"',
        "  fi",
        '  echo "$out" >&2',
        '  echo "${out##* }"',
        "}",
        "",
        "join_deps() {",
        "  local joined=''",
        '  for item in "$@"; do',
        '    if [[ -z "$item" ]]; then',
        "      continue",
        "    fi",
        '    if [[ -z "$joined" ]]; then',
        '      joined="$item"',
        "    else",
        '      joined="${joined}:$item"',
        "    fi",
        "  done",
        '  echo "$joined"',
        "}",
        "",
        'SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
        'prev_stage_deps=""',
        "",
    ]
    for stage_name, scripts in stage_scripts:
        if not scripts:
            continue
        lines.append(f'echo "[run.sbatch] stage={stage_name}"')
        if stage_name == "split":
            lines.append(
                f'jid_split_bams="$(submit_with_dep "$SCRIPT_DIR/{scripts[0].name}" "$prev_stage_deps")"'
            )
            if len(scripts) > 1:
                lines.append(
                    f'jid_split_sort="$(submit_with_dep "$SCRIPT_DIR/{scripts[1].name}" "$jid_split_bams")"'
                )
                lines.append('prev_stage_deps="$jid_split_sort"')
            else:
                lines.append('prev_stage_deps="$jid_split_bams"')
        else:
            job_vars: list[str] = []
            for index, script_path in enumerate(scripts):
                var_name = f"jid_{stage_name}_{index}".replace("-", "_")
                lines.append(
                    f'{var_name}="$(submit_with_dep "$SCRIPT_DIR/{script_path.name}" "$prev_stage_deps")"'
                )
                job_vars.append(var_name)
            deps_join = " ".join(f"${var_name}" for var_name in job_vars)
            lines.append(f'prev_stage_deps="$(join_deps {deps_join})"')
        lines.append("")
    lines.append('echo "[run.sbatch] done final_dep=${prev_stage_deps}"')
    driver_command = "\n".join(lines)
    slurm_args = argparse.Namespace(
        job_name=f"dbit_all_driver_{settings['sample_id']}",
        slurm_partition=settings["slurm_partition"],
        slurm_mem=settings["slurm_mem"],
        slurm_cpus_per_task=settings["slurm_cpus_per_task"],
        slurm_output=settings["slurm_output"].replace(
            "%x", f"dbit_all_driver_{settings['sample_id']}"
        ),
        slurm_error=settings["slurm_error"].replace(
            "%x", f"dbit_all_driver_{settings['sample_id']}"
        ),
        module_line="",
    )
    generate_slurm_script(driver_command, output_path, log_dir, slurm_args)


def submit_generated_stage_scripts(
    generated_scripts: list[Path],
    runner: str,
    stage: str,
) -> str | None:
    if runner == "local":
        for script_path in generated_scripts:
            submit_script(script_path, runner)
        return None

    if stage == "split":
        if len(generated_scripts) < 2:
            raise ValueError("missing split stage scripts for slurm submission")
        jid_split_bams = submit_slurm_script(generated_scripts[0])
        jid_split_sort = submit_slurm_script(generated_scripts[1], jid_split_bams)
        return jid_split_sort

    submitted_job_ids = [submit_slurm_script(script_path) for script_path in generated_scripts]
    return join_dependency_ids(submitted_job_ids)


def main() -> int:
    args = parse_args()
    settings = resolve_settings(args)
    sample_work = Path(settings["work_root"]) / settings["sample_id"]
    command_dir = sample_work / "commands"
    log_dir = sample_work / "logs"
    split_submit_helper_path: Path | None = None

    if settings["stage"] != "all":
        validate_inputs_for_stage(
            settings["stage"],
            settings,
            sample_work,
            skip_workdir_inputs=bool(args.skip_workdir_input_checks),
        )

    if settings["stage"] == "all":
        # Same checks as each per-stage subprocess (--skip-workdir-input-checks), but in
        # the parent so bad paths fail with a clear ValueError before any subprocess.
        for stage_name in STAGE_SEQUENCE:
            validate_inputs_for_stage(
                stage_name,
                settings,
                sample_work,
                skip_workdir_inputs=True,
            )
        passthrough_args = build_stage_passthrough_args(sys.argv[1:])
        stage_scripts: list[tuple[str, list[Path]]] = []
        for stage_name in STAGE_SEQUENCE:
            stage_argv = [
                sys.executable,
                __file__,
                *passthrough_args,
                "--runner",
                settings["runner"],
                "--stage",
                stage_name,
            ]
            if settings["dry_run"]:
                stage_argv.append("--dry-run")
            stage_argv.append("--skip-workdir-input-checks")
            completed = subprocess.run(
                stage_argv, check=False, capture_output=True, text=True
            )
            if completed.stdout:
                print(completed.stdout, end="")
            if completed.stderr:
                print(completed.stderr, end="", file=sys.stderr)
            if completed.returncode != 0:
                sys.exit(completed.returncode)
            stage_scripts.append((stage_name, parse_generated_paths(completed.stdout)))

        driver_path: Path
        if settings["runner"] == "local":
            driver_path = command_dir / "run.sh"
            print(f"[make_cmd] script={driver_path}")
            if not settings["dry_run"]:
                generate_local_driver_script(stage_scripts, driver_path)
        else:
            driver_path = command_dir / "run.sbatch"
            print(f"[make_cmd] script={driver_path}")
            if not settings["dry_run"]:
                driver_settings = apply_stage_slurm_settings(
                    settings, args, STAGE_SEQUENCE[-1]
                )
                generate_slurm_driver_script(
                    stage_scripts, driver_path, log_dir, driver_settings
                )

        if not settings["dry_run"] and driver_path.exists():
            print(f"[make_cmd] generated={driver_path}")
        if settings["submit"] and not settings["dry_run"]:
            if settings["runner"] == "slurm":
                subprocess.run(["bash", str(driver_path)], check=True)
                print("[make_cmd] submitted_driver=1")
                print("[make_cmd] submit_mode=client_side_sbatch_dag")
            else:
                submit_script(driver_path, settings["runner"])
                print("[make_cmd] submitted_driver=1")

        print("[make_cmd] stage=all helper generation complete")
        return 0

    generated_scripts: list[Path] = []
    if settings["stage"] == "fastp_split":
        base_name = "01_fastp_split"
        command_args = argparse.Namespace(
            r1=settings["r1"],
            r2=settings["r2"],
            fastp_threads=settings["fastp_threads"],
            number_of_split_parts=settings["number_of_split_parts"],
            fastp_bin=settings["fastp_bin"],
        )
        command = build_fastp_split_command(command_args, sample_work)
        if settings["runner"] == "local":
            script_path = command_dir / f"{base_name}.sh"
        else:
            script_path = command_dir / f"{base_name}.sbatch"

        print(f"[make_cmd] runner={settings['runner']}")
        print(f"[make_cmd] stage={settings['stage']}")
        print(f"[make_cmd] sample_id={settings['sample_id']}")
        print(f"[make_cmd] script={script_path}")
        print(f"[make_cmd] command={command}")

        if settings["dry_run"]:
            return 0

        if settings["runner"] == "local":
            generate_local_script(command, script_path)
        else:
            slurm_args = argparse.Namespace(
                job_name=f"dbit_fastp_split_{settings['sample_id']}",
                slurm_partition=settings["slurm_partition"],
                slurm_mem=settings["slurm_mem"],
                slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                slurm_output=settings["slurm_output"],
                slurm_error=settings["slurm_error"],
                module_line="",
            )
            generate_slurm_script(command, script_path, log_dir, slurm_args)
        generated_scripts.append(script_path)
    elif settings["stage"] == "demux_extract_bc":
        command_args = argparse.Namespace(
            barcode1_whitelist=settings["barcode1_whitelist"],
            barcode2_whitelist=settings["barcode2_whitelist"],
            linker_bc=settings["linker_bc"],
            insert_left=settings["insert_left"],
            linker_edit_distance=settings["linker_edit_distance"],
            barcode_hamming_distance=settings["barcode_hamming_distance"],
            gzip_level=settings["gzip_level"],
            spike_edit_distance=settings["spike_edit_distance"],
            require_c_all_t=settings["require_c_all_t"],
        )
        if settings["runner"] == "local":
            script_path = command_dir / "02_demux_extract_bc.sh"
            command = build_demux_local_batch_command(command_args, sample_work)
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={script_path}")
            print(f"[make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            generated_scripts.append(script_path)
        else:
            chunks = build_demux_chunks_from_config(
                sample_work, settings["number_of_split_parts"]
            )
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] chunk_count={len(chunks)}")
            for chunk, r1_path, r2_path, out_prefix in chunks:
                base_name = f"02_demux_extract_bc_{chunk}"
                script_path = command_dir / f"{base_name}.sbatch"
                command = build_demux_chunk_command(
                    command_args, r1_path, r2_path, out_prefix
                )
                chunk_output = settings["slurm_output"].replace(
                    "%x", f"dbit_demux_{settings['sample_id']}_{chunk}"
                )
                chunk_error = settings["slurm_error"].replace(
                    "%x", f"dbit_demux_{settings['sample_id']}_{chunk}"
                )
                print(f"[make_cmd] script={script_path}")
                print(f"[make_cmd] command={command}")
                if not settings["dry_run"]:
                    slurm_args = argparse.Namespace(
                        job_name=f"dbit_demux_{settings['sample_id']}_{chunk}",
                        slurm_partition=settings["slurm_partition"],
                        slurm_mem=settings["slurm_mem"],
                        slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                        slurm_output=chunk_output,
                        slurm_error=chunk_error,
                        module_line="",
                    )
                    generate_slurm_script(command, script_path, log_dir, slurm_args)
                generated_scripts.append(script_path)
    elif settings["stage"] == "align":
        command_args = argparse.Namespace(
            bwa_index=settings["bwa_index"],
            bwa_threads=settings["bwa_threads"],
            bwa_bin=settings["bwa_bin"],
            sinto_bin=settings["sinto_bin"],
            samtools_bin=settings["samtools_bin"],
            spike_in_index=settings["spike_in_index"],
        )
        if settings["runner"] == "local":
            script_path = command_dir / "03_align.sh"
            command = build_align_local_batch_command(command_args, sample_work)
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={script_path}")
            print(f"[make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            generated_scripts.append(script_path)
        else:
            chunks = build_align_chunks_from_config(
                sample_work, settings["number_of_split_parts"]
            )
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] chunk_count={len(chunks)}")
            for chunk, _, _ in chunks:
                base_name = f"03_align_{chunk}"
                script_path = command_dir / f"{base_name}.sbatch"
                command = build_align_chunk_command(command_args, sample_work, chunk)
                chunk_output = settings["slurm_output"].replace(
                    "%x", f"dbit_align_{settings['sample_id']}_{chunk}"
                )
                chunk_error = settings["slurm_error"].replace(
                    "%x", f"dbit_align_{settings['sample_id']}_{chunk}"
                )
                print(f"[make_cmd] script={script_path}")
                print(f"[make_cmd] command={command}")
                if not settings["dry_run"]:
                    slurm_args = argparse.Namespace(
                        job_name=f"dbit_align_{settings['sample_id']}_{chunk}",
                        slurm_partition=settings["slurm_partition"],
                        slurm_mem=settings["slurm_mem"],
                        slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                        slurm_output=chunk_output,
                        slurm_error=chunk_error,
                        module_line="",
                    )
                    generate_slurm_script(command, script_path, log_dir, slurm_args)
                generated_scripts.append(script_path)
    elif settings["stage"] == "pool":
        command_args = argparse.Namespace(
            samtools_bin=settings["samtools_bin"],
            samtools_threads=settings["samtools_threads"],
            host_sort_mem=settings["host_sort_mem"],
            spike_in_index=settings["spike_in_index"],
        )
        if settings["runner"] == "local":
            script_path = command_dir / "04_pool.sh"
            command = build_pool_command(command_args, sample_work, "all")
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={script_path}")
            print(f"[make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            generated_scripts.append(script_path)
        else:
            host_script_path = command_dir / "04_pool_host.sbatch"
            spike_script_path = command_dir / "04_pool_spike.sbatch"
            host_command = build_pool_command(command_args, sample_work, "host")
            spike_command = build_pool_command(command_args, sample_work, "spike")
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={spike_script_path}")
            print(f"[make_cmd] command={spike_command}")
            print(f"[make_cmd] script={host_script_path}")
            print(f"[make_cmd] command={host_command}")
            if not settings["dry_run"]:
                spike_output = settings["slurm_output"].replace(
                    "%x", f"dbit_pool_spike_{settings['sample_id']}"
                )
                spike_error = settings["slurm_error"].replace(
                    "%x", f"dbit_pool_spike_{settings['sample_id']}"
                )
                host_output = settings["slurm_output"].replace(
                    "%x", f"dbit_pool_host_{settings['sample_id']}"
                )
                host_error = settings["slurm_error"].replace(
                    "%x", f"dbit_pool_host_{settings['sample_id']}"
                )
                spike_slurm_args = argparse.Namespace(
                    job_name=f"dbit_pool_spike_{settings['sample_id']}",
                    slurm_partition=settings["slurm_partition"],
                    slurm_mem=settings["slurm_mem"],
                    slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                    slurm_output=spike_output,
                    slurm_error=spike_error,
                    module_line="",
                )
                host_slurm_args = argparse.Namespace(
                    job_name=f"dbit_pool_host_{settings['sample_id']}",
                    slurm_partition=settings["slurm_partition"],
                    slurm_mem=settings["slurm_mem"],
                    slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                    slurm_output=host_output,
                    slurm_error=host_error,
                    module_line="",
                )
                generate_slurm_script(
                    spike_command, spike_script_path, log_dir, spike_slurm_args
                )
                generate_slurm_script(host_command, host_script_path, log_dir, host_slurm_args)
            generated_scripts.append(spike_script_path)
            generated_scripts.append(host_script_path)
    elif settings["stage"] == "split":
        command_args = argparse.Namespace(
            split_barcodes=settings["split_barcodes"],
            split_cb_tag=settings["split_cb_tag"],
            split_threads_read=settings["split_threads_read"],
            split_threads_write=settings["split_threads_write"],
            split_smoke=settings["split_smoke"],
            samtools_bin=settings["samtools_bin"],
            split_sort_jobs=settings["split_sort_jobs"],
        )
        if settings["runner"] == "local":
            script_path = command_dir / "05_split.sh"
            split_command = build_split_command(command_args, sample_work)
            sort_command = build_split_sort_command(command_args, sample_work)
            command = f"{split_command}\n{sort_command}"
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={script_path}")
            print(f"[make_cmd] command={split_command}")
            print(f"[make_cmd] command={sort_command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            generated_scripts.append(script_path)
        else:
            split_script_path = command_dir / "05_split_bams.sbatch"
            sort_script_path = command_dir / "05_split_sort.sbatch"
            split_submit_helper_path = command_dir / "05_split_submit.sh"
            split_command = build_split_command(command_args, sample_work)
            sort_command = build_split_sort_command(command_args, sample_work)
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={split_script_path}")
            print(f"[make_cmd] command={split_command}")
            print(f"[make_cmd] script={sort_script_path}")
            print(f"[make_cmd] command={sort_command}")
            print(f"[make_cmd] helper={split_submit_helper_path}")
            if not settings["dry_run"]:
                split_slurm_args = argparse.Namespace(
                    job_name=f"dbit_split_bams_{settings['sample_id']}",
                    slurm_partition=settings["split_bams_slurm_partition"],
                    slurm_mem=settings["split_bams_slurm_mem"],
                    slurm_cpus_per_task=settings["split_bams_slurm_cpus_per_task"],
                    slurm_output=settings["split_bams_slurm_output"].replace(
                        "%x", f"dbit_split_bams_{settings['sample_id']}"
                    ),
                    slurm_error=settings["split_bams_slurm_error"].replace(
                        "%x", f"dbit_split_bams_{settings['sample_id']}"
                    ),
                    module_line="",
                )
                sort_slurm_args = argparse.Namespace(
                    job_name=f"dbit_split_sort_{settings['sample_id']}",
                    slurm_partition=settings["split_sort_slurm_partition"],
                    slurm_mem=settings["split_sort_slurm_mem"],
                    slurm_cpus_per_task=settings["split_sort_slurm_cpus_per_task"],
                    slurm_output=settings["split_sort_slurm_output"].replace(
                        "%x", f"dbit_split_sort_{settings['sample_id']}"
                    ),
                    slurm_error=settings["split_sort_slurm_error"].replace(
                        "%x", f"dbit_split_sort_{settings['sample_id']}"
                    ),
                    module_line="",
                )
                generate_slurm_script(
                    split_command, split_script_path, log_dir, split_slurm_args
                )
                generate_slurm_script(
                    sort_command, sort_script_path, log_dir, sort_slurm_args
                )
                generate_split_slurm_submit_script(
                    split_script_path, sort_script_path, split_submit_helper_path
                )
            generated_scripts.append(split_script_path)
            generated_scripts.append(sort_script_path)
    elif settings["stage"] == "mbias":
        command_args = argparse.Namespace(
            samtools_bin=settings["samtools_bin"],
            samtools_threads=settings["samtools_threads"],
            mbias_host_subsample_fraction=settings["mbias_host_subsample_fraction"],
            mbias_max_cycle=settings["mbias_max_cycle"],
            mbias_min_mapping_quality=settings["mbias_min_mapping_quality"],
            mbias_script=settings["mbias_script"],
            call_reference_file=settings["call_reference_file"],
            call_chromosomes=settings["call_chromosomes"],
            spike_in_index=settings["spike_in_index"],
        )
        if settings["runner"] == "local":
            script_path = command_dir / "06_mbias.sh"
            command = build_mbias_command(
                command_args,
                sample_work,
                settings["mbias_mode"],
            )
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={script_path}")
            print(f"[make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            generated_scripts.append(script_path)
        else:
            mbias_mode = settings["mbias_mode"]
            include_host = mbias_mode in ("all", "host")
            include_spike = mbias_mode in ("all", "spike")
            host_script_path = command_dir / "06_mbias_host.sbatch"
            host_command = build_mbias_command(command_args, sample_work, "host")
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            if include_host:
                print(f"[make_cmd] script={host_script_path}")
                print(f"[make_cmd] command={host_command}")
                if not settings["dry_run"]:
                    host_output = settings["mbias_host_slurm_output"].replace(
                        "%x", f"dbit_mbias_host_{settings['sample_id']}"
                    )
                    host_error = settings["mbias_host_slurm_error"].replace(
                        "%x", f"dbit_mbias_host_{settings['sample_id']}"
                    )
                    host_slurm_args = argparse.Namespace(
                        job_name=f"dbit_mbias_host_{settings['sample_id']}",
                        slurm_partition=settings["mbias_host_slurm_partition"],
                        slurm_mem=settings["mbias_host_slurm_mem"],
                        slurm_cpus_per_task=settings["mbias_host_slurm_cpus_per_task"],
                        slurm_output=host_output,
                        slurm_error=host_error,
                        module_line="",
                    )
                    generate_slurm_script(host_command, host_script_path, log_dir, host_slurm_args)
                generated_scripts.append(host_script_path)
            if include_spike:
                spike_names = parse_spike_names(settings["spike_in_index"])
                for spike_name in spike_names:
                    spike_script_path = command_dir / f"06_mbias_spike_{spike_name}.sbatch"
                    spike_command = build_mbias_command(
                        command_args,
                        sample_work,
                        "spike",
                        spike_name=spike_name,
                    )
                    print(f"[make_cmd] script={spike_script_path}")
                    print(f"[make_cmd] command={spike_command}")
                    if not settings["dry_run"]:
                        spike_output = settings["mbias_spike_slurm_output"].replace(
                            "%x", f"dbit_mbias_spike_{settings['sample_id']}_{spike_name}"
                        )
                        spike_error = settings["mbias_spike_slurm_error"].replace(
                            "%x", f"dbit_mbias_spike_{settings['sample_id']}_{spike_name}"
                        )
                        spike_slurm_args = argparse.Namespace(
                            job_name=f"dbit_mbias_spike_{settings['sample_id']}_{spike_name}",
                            slurm_partition=settings["mbias_spike_slurm_partition"],
                            slurm_mem=settings["mbias_spike_slurm_mem"],
                            slurm_cpus_per_task=settings["mbias_spike_slurm_cpus_per_task"],
                            slurm_output=spike_output,
                            slurm_error=spike_error,
                            module_line="",
                        )
                        generate_slurm_script(
                            spike_command, spike_script_path, log_dir, spike_slurm_args
                        )
                    generated_scripts.append(spike_script_path)
    elif settings["stage"] == "call":
        command_args = argparse.Namespace(
            call_reference_file=settings["call_reference_file"],
            call_chromosomes=settings["call_chromosomes"],
            call_mito_chromosomes=settings["call_mito_chromosomes"],
            call_jobs=settings["call_jobs"],
            call_min_base_quality=settings["call_min_base_quality"],
            call_min_mapping_quality=settings["call_min_mapping_quality"],
            call_sample_size=settings["call_sample_size"],
            call_max_depth=settings["call_max_depth"],
            call_batch_size=settings["call_batch_size"],
            call_context_mode=settings["call_context_mode"],
            call_trim_host=settings["call_trim_host"],
            call_trim_spike=settings["call_trim_spike"],
            call_caller_script=settings["call_caller_script"],
            call_ch_caller_script=settings["call_ch_caller_script"],
            samtools_bin=settings["samtools_bin"],
            samtools_threads=settings["samtools_threads"],
            mbias_host_subsample_fraction=settings["mbias_host_subsample_fraction"],
            spike_in_index=settings["spike_in_index"],
        )
        if settings["runner"] == "local":
            script_path = command_dir / "07_call.sh"
            call_mode = settings["call_mode"]
            trims_differ = (
                settings["call_trim_host"] != settings["call_trim_spike"]
            )
            if call_mode == "all" and trims_differ:
                sub_commands: list[str] = [
                    build_call_command(command_args, sample_work, "host")
                ]
                for spike_name in parse_spike_names(settings["spike_in_index"]):
                    sub_commands.append(
                        build_call_command(
                            command_args,
                            sample_work,
                            "spike",
                            spike_name=spike_name,
                        )
                    )
                command = "\n\n".join(sub_commands)
            else:
                command = build_call_command(command_args, sample_work, call_mode)
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={script_path}")
            print(f"[make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            generated_scripts.append(script_path)
        else:
            call_mode = settings["call_mode"]
            include_host = call_mode in ("all", "host")
            include_spike = call_mode in ("all", "spike")
            host_script_path = command_dir / "07_call_host.sbatch"
            host_command = build_call_command(command_args, sample_work, "host")
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            if include_host:
                print("[make_cmd] host_spot_count=runtime_discovery")
                print(f"[make_cmd] script={host_script_path}")
                print(f"[make_cmd] command={host_command}")
                if not settings["dry_run"]:
                    host_output = settings["call_host_slurm_output"].replace(
                        "%x", f"dbit_call_host_{settings['sample_id']}"
                    )
                    host_error = settings["call_host_slurm_error"].replace(
                        "%x", f"dbit_call_host_{settings['sample_id']}"
                    )
                    host_slurm_args = argparse.Namespace(
                        job_name=f"dbit_call_host_{settings['sample_id']}",
                        slurm_partition=settings["call_host_slurm_partition"],
                        slurm_mem=settings["call_host_slurm_mem"],
                        slurm_cpus_per_task=settings["call_host_slurm_cpus_per_task"],
                        slurm_output=host_output,
                        slurm_error=host_error,
                        module_line="",
                    )
                    generate_slurm_script(host_command, host_script_path, log_dir, host_slurm_args)
                generated_scripts.append(host_script_path)

            if include_spike:
                spike_names = parse_spike_names(settings["spike_in_index"])
                for spike_name in spike_names:
                    spike_script_path = command_dir / f"07_call_spike_{spike_name}.sbatch"
                    spike_command = build_call_command(
                        command_args,
                        sample_work,
                        "spike",
                        spike_name=spike_name,
                    )
                    print(f"[make_cmd] script={spike_script_path}")
                    print(f"[make_cmd] command={spike_command}")
                    if not settings["dry_run"]:
                        spike_output = settings["call_spike_slurm_output"].replace(
                            "%x", f"dbit_call_spike_{settings['sample_id']}_{spike_name}"
                        )
                        spike_error = settings["call_spike_slurm_error"].replace(
                            "%x", f"dbit_call_spike_{settings['sample_id']}_{spike_name}"
                        )
                        spike_slurm_args = argparse.Namespace(
                            job_name=f"dbit_call_spike_{settings['sample_id']}_{spike_name}",
                            slurm_partition=settings["call_spike_slurm_partition"],
                            slurm_mem=settings["call_spike_slurm_mem"],
                            slurm_cpus_per_task=settings["call_spike_slurm_cpus_per_task"],
                            slurm_output=spike_output,
                            slurm_error=spike_error,
                            module_line="",
                        )
                        generate_slurm_script(
                            spike_command, spike_script_path, log_dir, spike_slurm_args
                        )
                    generated_scripts.append(spike_script_path)
    elif settings["stage"] == "saturation":
        command_args = argparse.Namespace(
            saturation_script=settings["saturation_script"],
            saturation_reads_threshold=settings["saturation_reads_threshold"],
        )
        command = build_saturation_command(command_args, sample_work)
        if settings["runner"] == "local":
            script_path = command_dir / "08_saturation.sh"
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={script_path}")
            print(f"[make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            generated_scripts.append(script_path)
        else:
            script_path = command_dir / "08_saturation.sbatch"
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={script_path}")
            print(f"[make_cmd] command={command}")
            if not settings["dry_run"]:
                slurm_args = argparse.Namespace(
                    job_name=f"dbit_saturation_{settings['sample_id']}",
                    slurm_partition=settings["slurm_partition"],
                    slurm_mem=settings["slurm_mem"],
                    slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                    slurm_output=settings["slurm_output"].replace(
                        "%x", f"dbit_saturation_{settings['sample_id']}"
                    ),
                    slurm_error=settings["slurm_error"].replace(
                        "%x", f"dbit_saturation_{settings['sample_id']}"
                    ),
                    module_line="",
                )
                generate_slurm_script(command, script_path, log_dir, slurm_args)
            generated_scripts.append(script_path)
    elif settings["stage"] == "summary":
        command_args = argparse.Namespace(
            summary_script=settings["summary_script"],
            spike_in_index=settings["spike_in_index"],
            call_mito_chromosomes=settings["call_mito_chromosomes"],
        )
        command = build_summary_command(command_args, sample_work)
        if settings["runner"] == "local":
            script_path = command_dir / "09_summary.sh"
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={script_path}")
            print(f"[make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            generated_scripts.append(script_path)
        else:
            script_path = command_dir / "09_summary.sbatch"
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={script_path}")
            print(f"[make_cmd] command={command}")
            if not settings["dry_run"]:
                slurm_args = argparse.Namespace(
                    job_name=f"dbit_summary_{settings['sample_id']}",
                    slurm_partition=settings["slurm_partition"],
                    slurm_mem=settings["slurm_mem"],
                    slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                    slurm_output=settings["slurm_output"].replace(
                        "%x", f"dbit_summary_{settings['sample_id']}"
                    ),
                    slurm_error=settings["slurm_error"].replace(
                        "%x", f"dbit_summary_{settings['sample_id']}"
                    ),
                    module_line="",
                )
                generate_slurm_script(command, script_path, log_dir, slurm_args)
            generated_scripts.append(script_path)
    elif settings["stage"] == "aggregate":
        command_args = argparse.Namespace(
            aggregate_script=settings["aggregate_script"],
            sort_mem=settings["aggregate_sort_mem"],
        )
        command = build_aggregate_command(command_args, sample_work)
        if settings["runner"] == "local":
            script_path = command_dir / "10_aggregate.sh"
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={script_path}")
            print(f"[make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            generated_scripts.append(script_path)
        else:
            script_path = command_dir / "10_aggregate.sbatch"
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={script_path}")
            print(f"[make_cmd] command={command}")
            if not settings["dry_run"]:
                slurm_args = argparse.Namespace(
                    job_name=f"dbit_aggregate_{settings['sample_id']}",
                    slurm_partition=settings["slurm_partition"],
                    slurm_mem=settings["slurm_mem"],
                    slurm_cpus_per_task=settings["slurm_cpus_per_task"],
                    slurm_output=settings["slurm_output"].replace(
                        "%x", f"dbit_aggregate_{settings['sample_id']}"
                    ),
                    slurm_error=settings["slurm_error"].replace(
                        "%x", f"dbit_aggregate_{settings['sample_id']}"
                    ),
                    module_line="",
                )
                generate_slurm_script(command, script_path, log_dir, slurm_args)
            generated_scripts.append(script_path)
    elif settings["stage"] in TAPS_METHSCAN_STAGE_MAP:
        step, script_base, slug = TAPS_METHSCAN_STAGE_MAP[settings["stage"]]
        ext = ".sh" if settings["runner"] == "local" else ".sbatch"
        job_name = f"{slug}_{settings['sample_id']}"
        emit_methscan_stage_command(
            settings=settings,
            sample_work=sample_work,
            command_dir=command_dir,
            log_dir=log_dir,
            step=step,
            script_filename=f"{script_base}{ext}",
            job_name=job_name,
            generated_scripts=generated_scripts,
        )
        if settings["dry_run"]:
            return 0
    else:
        raise ValueError(f"unsupported stage: {settings['stage']}")

    for script_path in generated_scripts:
        print(f"[make_cmd] generated={script_path}")
    if split_submit_helper_path and not settings["dry_run"] and split_submit_helper_path.exists():
        print(f"[make_cmd] helper_generated={split_submit_helper_path}")

    if settings["submit"]:
        final_dependency = submit_generated_stage_scripts(
            generated_scripts,
            settings["runner"],
            settings["stage"],
        )
        if final_dependency:
            print(f"[make_cmd] final_slurm_dependency={final_dependency}")
        print(f"[make_cmd] submitted_count={len(generated_scripts)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
