#!/usr/bin/env python3
"""Generate and optionally submit EMSeq workflow commands.

This is an EMSeq-only entrypoint that supports:

- `fastp_split`
- `demux_extract_bc`
- `align`
- `pool`
- `split`
- `mbias`
- `call`
- `saturation`
- `summary`
- `aggregate` (experimental; optional, not part of `--stage all`; reuses `scripts/aggregate.py`)
- `methscan_prepare` (experimental; optional, not part of `--stage all`; host `*.CG.cov` via `scripts/methscan_prepare.py`)
- `all` (generates `commands/run.sh` or `commands/run.sbatch` to run the full pipeline in order)

It intentionally keeps EMSeq orchestration separate from the TAPS workflow
while reusing only the stage implementations whose contracts still match the
EMSeq library layout.
"""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
from pathlib import Path


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
# Stages that may appear as top-level keys under workflow `slurm` (nested mode).
SLURM_NEST_STAGE_KEYS = frozenset(STAGE_SEQUENCE) | {"aggregate", "methscan_prepare"}
STAGE_CHOICES = [*STAGE_SEQUENCE, "aggregate", "methscan_prepare", "all"]
STAGE_REQUIRED_FIELDS = {
    "fastp_split": ["r1", "r2", "number_of_split_parts"],
    "demux_extract_bc": [
        "barcode1_whitelist",
        "barcode2_whitelist",
        "number_of_split_parts",
    ],
    "align": [
        "number_of_split_parts",
        "biscuit_reference",
    ],
    "pool": [],
    "split": ["split_barcodes"],
    "mbias": [],
    "call": ["call_reference_file", "call_jobs"],
    "saturation": [],
    "summary": [],
    "aggregate": [],
    "methscan_prepare": [],
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
        description="Generate executable command scripts for EMSeq workflow stages."
    )
    parser.add_argument(
        "--workflow-config",
        help="JSON config path for EMSeq workflow/sample settings.",
    )
    parser.add_argument(
        "--runner",
        choices=["local", "slurm"],
        help="Command target: local shell or slurm sbatch.",
    )
    parser.add_argument(
        "--stage",
        choices=STAGE_CHOICES,
        help=(
            "Workflow stage to generate command script for. "
            "Use `all` for a driver that runs every stage in order. "
            "Default: fastp_split."
        ),
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
        "--biscuit-reference",
        help="Reference FASTA passed to biscuit align (EMSeq align stage).",
    )
    parser.add_argument(
        "--biscuit-threads",
        type=int,
        help="Thread count for biscuit align -@. Default: 32 in EMSeq config.",
    )
    parser.add_argument(
        "--biscuit-batch-size",
        type=int,
        help="Batch size for biscuit align -b. Default: 1.",
    )
    parser.add_argument(
        "--biscuit-bin",
        help=(
            "biscuit executable path or command name. "
            "Default: biscuit from current Python env if available, else biscuit."
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
    # Pool stage.
    parser.add_argument(
        "--samtools-threads",
        type=int,
        help="Thread count for samtools sort in pool stage. Default: 4.",
    )
    parser.add_argument(
        "--host-sort-mem",
        help="Memory per thread for host samtools sort -m in pool stage. Default: 16G.",
    )
    # Split stage.
    parser.add_argument(
        "--split-barcodes",
        help="Barcode whitelist TSV for split stage. Default: barcode1_whitelist.",
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
        help="Enable smoke mode for split stage (emit up to 16 spot BAMs).",
    )
    parser.add_argument(
        "--split-sort-jobs",
        type=int,
        help="Parallel job count for bam_sort_parallel in split stage. Default: 8.",
    )
    # Call stage.
    parser.add_argument(
        "--call-reference-file",
        help="Host reference FASTA for EMSeq call stage (biscuit pileup/mergecg).",
    )
    parser.add_argument(
        "--call-jobs",
        type=int,
        help="Maximum per-stage pileup jobs for EMSeq call. Default: 8.",
    )
    parser.add_argument(
        "--call-mode",
        choices=["all", "host", "spike"],
        help="Call mode: host+spike, host only, or spike only. Default: all.",
    )
    parser.add_argument(
        "--call-host-threads",
        type=int,
        help="biscuit pileup/bgzip threads in host call. Default: biscuit_threads.",
    )
    parser.add_argument(
        "--call-spike-threads",
        type=int,
        help="biscuit pileup/bgzip threads in spike call. Default: biscuit_threads.",
    )
    parser.add_argument(
        "--call-bgzip-bin",
        help="bgzip executable path or command name. Default: bgzip.",
    )
    parser.add_argument(
        "--call-tabix-bin",
        help="tabix executable path or command name. Default: tabix.",
    )
    parser.add_argument(
        "--call-host-subsample-fraction",
        type=float,
        help="Host subsample fraction for host_mito BAM fallback. Default: 0.1.",
    )
    parser.add_argument(
        "--call-host-subsample-seed",
        type=int,
        help="Host subsample seed for host_mito BAM fallback. Default: scripts.host_subsample_bam.HOST_SUBSAMPLE_SEED.",
    )
    parser.add_argument(
        "--call-mito-chromosomes",
        help=(
            "Comma-separated host contigs treated as mitochondrial in EMSeq call. "
            "Default: chrM."
        ),
    )
    parser.add_argument(
        "--call-left-trimming",
        type=int,
        default=None,
        metavar="N",
        help=(
            "If set, biscuit pileup -5 N (read 5' end distance). "
            "If omitted, use workflow default or biscuit built-in default."
        ),
    )
    parser.add_argument(
        "--call-right-trimming",
        type=int,
        default=None,
        metavar="N",
        help=(
            "If set, biscuit pileup -3 N (read 3' end distance). "
            "If omitted, use workflow default or biscuit built-in default."
        ),
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
    # M-bias stage (scripts-emseq/mbias.py; runs after split, before call).
    parser.add_argument(
        "--mbias-script",
        help=(
            "Path to EMSeq mbias script. "
            "Default: scripts-emseq/mbias.py."
        ),
    )
    parser.add_argument(
        "--mbias-mode",
        choices=["all", "host", "spike"],
        help="mbias mode: host+spike, host only, or spike only. Default: spike.",
    )
    parser.add_argument(
        "--mbias-host-subsample-fraction",
        type=float,
        help="Host subsampling fraction for mbias stage. Default: 0.1.",
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
    # Saturation stage (reuses scripts/saturation.py; runs after call).
    parser.add_argument(
        "--saturation-script",
        help=(
            "Path to saturation script. "
            "Default: scripts/saturation.py."
        ),
    )
    parser.add_argument(
        "--saturation-reads-threshold",
        type=float,
        help=(
            "HQ spot threshold (reads) for scripts/saturation.py. "
            "Default: 1e6."
        ),
    )
    # Summary stage (reuses scripts/summary.py; runs after saturation).
    parser.add_argument(
        "--summary-script",
        help=(
            "Path to summary script. "
            "Default: scripts/summary.py."
        ),
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
        "--methscan-prepare-script",
        help=(
            "Path to methscan_prepare wrapper. Default: scripts/methscan_prepare.py."
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
        "--methscan-prepare-chunksize",
        type=int,
        default=None,
        metavar="N",
        help=(
            "Passed to scripts/methscan_prepare.py --chunksize (methscan prepare). "
            "Default: 10000000 (10 Mbp)."
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
    script_path = Path("scripts/extract_bc.py")
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
        ]
    )


def build_demux_local_batch_command(
    args: argparse.Namespace, sample_work: Path
) -> str:
    chunk_dir = sample_work / "shard_fastq"
    demux_dir = sample_work / "demux"
    chunk_dir_q = shlex.quote(str(chunk_dir))
    demux_dir_q = shlex.quote(str(demux_dir))
    py = quoted([sys.executable, "scripts/extract_bc.py"])
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
        f"--gzip-level {int(args.gzip_level)}\n"
        "done\n"
        "\n"
        'echo "[demux] done"'
    )


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


def parse_spike_names(spike_in_index: list[str]) -> list[str]:
    """Extract unique spike-in names from NAME=INDEX items."""
    names: list[str] = []
    for item in spike_in_index:
        name, _, _ = item.partition("=")
        name = name.strip()
        if name and name not in names:
            names.append(name)
    return names


def build_pool_command(args: argparse.Namespace, sample_work: Path, mode: str) -> str:
    """Build scripts/pool.py invocation."""
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
        str(args.host_sort_mem),
        "--mode",
        mode,
    ]
    for spike_name in parse_spike_names(args.spike_in_index):
        command.extend(["--spike-in-name", spike_name])
    return quoted(command)


def build_split_command(args: argparse.Namespace, sample_work: Path) -> str:
    """Build scripts/split_bams.py invocation."""
    script_path = Path("scripts/split_bams.py")
    command = [
        sys.executable,
        str(script_path),
        "--in-bam",
        str(sample_work / "pooled" / "pooled.byCB.bam"),
        "--barcodes",
        str(args.split_barcodes),
        "--out-dir",
        str(sample_work / "split_bams"),
        "--cb-tag",
        str(args.split_cb_tag),
        "--threads-read",
        str(args.split_threads_read),
        "--threads-write",
        str(args.split_threads_write),
    ]
    if args.split_smoke:
        command.append("--smoke")
    return quoted(command)


def build_split_sort_command(args: argparse.Namespace, sample_work: Path) -> str:
    """Build scripts/bam_sort_parallel.py invocation."""
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


def build_mbias_command(
    args: argparse.Namespace,
    sample_work: Path,
    mode: str,
    spike_name: str | None = None,
) -> str:
    """Build scripts-emseq/mbias.py invocation."""
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
    if mode in ("all", "spike"):
        for item in args.spike_in_index:
            command.extend(["--spike-reference", item])
    if spike_name:
        command.extend(["--spike-in-name", spike_name])
    return quoted(command)


def build_call_command(
    args: argparse.Namespace,
    sample_work: Path,
    mode: str,
    spike_name: str | None = None,
    spike_references: list[str] | None = None,
) -> str:
    script_path = Path("scripts-emseq/call.py")
    command: list[str] = [
        sys.executable,
        str(script_path),
        "--work-path",
        str(sample_work),
        "--mode",
        mode,
        "--reference-file",
        args.call_reference_file,
        "--jobs",
        str(args.call_jobs),
        "--host-threads",
        str(args.call_host_threads),
        "--spike-threads",
        str(args.call_spike_threads),
        "--biscuit-bin",
        args.biscuit_bin,
        "--bgzip-bin",
        args.call_bgzip_bin,
        "--tabix-bin",
        args.call_tabix_bin,
        "--mito-chromosomes",
        args.call_mito_chromosomes,
    ]
    if args.call_left_trimming is not None:
        command.extend(["--call-left-trimming", str(args.call_left_trimming)])
    if args.call_right_trimming is not None:
        command.extend(["--call-right-trimming", str(args.call_right_trimming)])
    if spike_references:
        for item in spike_references:
            command.extend(["--spike-reference", item])
    if spike_name:
        command.extend(["--spike-in-name", spike_name])
    return quoted(command)


def build_saturation_command(args: argparse.Namespace, sample_work: Path) -> str:
    """Build scripts/saturation.py invocation."""
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


def build_summary_command(args: argparse.Namespace, sample_work: Path) -> str:
    """Build scripts/summary.py invocation (spike names from spike_in_index)."""
    script_path = Path(args.summary_script)
    command = [
        sys.executable,
        str(script_path),
        "--work-path",
        str(sample_work),
    ]
    for spike_name in parse_spike_names(args.spike_in_index):
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


def build_methscan_prepare_command(args: argparse.Namespace, sample_work: Path) -> str:
    script_path = Path(args.methscan_prepare_script)
    command = [
        sys.executable,
        str(script_path),
        "--work-path",
        str(sample_work),
        "--chunksize",
        str(args.methscan_prepare_chunksize),
    ]
    manifest = getattr(args, "methscan_pixi_manifest", None)
    if manifest:
        command.extend(["--pixi-manifest", str(manifest)])
    return quoted(command)


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


def select_stage_slurm_cfg(slurm_cfg_raw: dict, stage: str) -> dict:
    if any(key in slurm_cfg_raw for key in SLURM_NEST_STAGE_KEYS):
        stage_slurm_cfg = slurm_cfg_raw.get(stage, {})
    else:
        stage_slurm_cfg = slurm_cfg_raw
    if stage_slurm_cfg is None:
        stage_slurm_cfg = {}
    if not isinstance(stage_slurm_cfg, dict):
        raise ValueError("selected slurm config must be an object")
    return stage_slurm_cfg


def resolve_step_slurm_cfg(
    stage_slurm_cfg: dict,
    step_name: str,
    step_keys: set[str],
) -> dict:
    """Resolve nested slurm step configuration (e.g. slurm.split.split_bams)."""
    if not isinstance(stage_slurm_cfg, dict):
        raise ValueError("selected slurm config must be an object")
    if not any(key in stage_slurm_cfg for key in step_keys):
        return stage_slurm_cfg
    step_cfg = stage_slurm_cfg.get(step_name, {})
    if step_cfg is None:
        step_cfg = {}
    if not isinstance(step_cfg, dict):
        raise ValueError(f"slurm step config '{step_name}' must be an object")
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
    if stage not in STAGE_CHOICES:
        raise ValueError(f"unsupported stage for EMSeq entry: {stage}")

    if stage == "all":
        slurm_stage_key = STAGE_SEQUENCE[0]
    elif stage in STAGE_SEQUENCE or stage in ("aggregate", "methscan_prepare"):
        slurm_stage_key = stage
    else:
        slurm_stage_key = STAGE_SEQUENCE[0]
    stage_slurm_cfg = select_stage_slurm_cfg(slurm_cfg_raw, slurm_stage_key)

    split_bams_slurm_cfg = resolve_step_slurm_cfg(
        stage_slurm_cfg, "split_bams", {"split_bams", "sort"}
    )
    split_sort_slurm_cfg = resolve_step_slurm_cfg(
        stage_slurm_cfg, "sort", {"split_bams", "sort"}
    )

    mbias_stage_slurm_cfg = select_stage_slurm_cfg(slurm_cfg_raw, "mbias")
    mbias_host_slurm_cfg = resolve_step_slurm_cfg(
        mbias_stage_slurm_cfg, "host", {"host", "spike"}
    )
    mbias_spike_slurm_cfg = resolve_step_slurm_cfg(
        mbias_stage_slurm_cfg, "spike", {"host", "spike"}
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
        "biscuit_reference": pick(args.biscuit_reference, cfg.get("biscuit_reference")),
        "biscuit_threads": pick(args.biscuit_threads, cfg.get("biscuit_threads")),
        "biscuit_batch_size": pick(
            args.biscuit_batch_size, cfg.get("biscuit_batch_size")
        ),
        "biscuit_bin": pick(args.biscuit_bin, cfg.get("biscuit_bin")),
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
        # Call stage.
        "call_reference_file": pick(
            args.call_reference_file,
            cfg.get("call_reference_file", cfg.get("biscuit_reference")),
        ),
        "call_jobs": pick(args.call_jobs, cfg.get("call_jobs")),
        "call_mode": pick(args.call_mode, cfg.get("call_mode")),
        "call_host_threads": pick(
            args.call_host_threads, cfg.get("call_host_threads")
        ),
        "call_spike_threads": pick(
            args.call_spike_threads, cfg.get("call_spike_threads")
        ),
        "call_bgzip_bin": pick(args.call_bgzip_bin, cfg.get("call_bgzip_bin")),
        "call_tabix_bin": pick(args.call_tabix_bin, cfg.get("call_tabix_bin")),
        "host_subsample_fraction": pick(
            args.call_host_subsample_fraction, cfg.get("call_host_subsample_fraction")
        ),
        "host_subsample_seed": pick(
            args.call_host_subsample_seed, cfg.get("call_host_subsample_seed")
        ),
        "call_mito_chromosomes": pick(
            args.call_mito_chromosomes, cfg.get("call_mito_chromosomes")
        ),
        "call_left_trimming": pick(
            args.call_left_trimming, cfg.get("call_left_trimming")
        ),
        "call_right_trimming": pick(
            args.call_right_trimming, cfg.get("call_right_trimming")
        ),
        "spike_in_index": pick(args.spike_in_index, cfg.get("spike_in_index")),
        # M-bias stage.
        "mbias_script": pick(args.mbias_script, cfg.get("mbias_script")),
        "mbias_mode": pick(args.mbias_mode, cfg.get("mbias_mode")),
        "mbias_host_subsample_fraction": pick(
            args.mbias_host_subsample_fraction,
            cfg.get("mbias_host_subsample_fraction"),
        ),
        "mbias_max_cycle": pick(args.mbias_max_cycle, cfg.get("mbias_max_cycle")),
        "mbias_min_mapping_quality": pick(
            args.mbias_min_mapping_quality, cfg.get("mbias_min_mapping_quality")
        ),
        # Saturation stage.
        "saturation_script": pick(args.saturation_script, cfg.get("saturation_script")),
        "saturation_reads_threshold": pick(
            args.saturation_reads_threshold,
            cfg.get("saturation_reads_threshold"),
        ),
        # Summary stage.
        "summary_script": pick(args.summary_script, cfg.get("summary_script")),
        "aggregate_script": pick(args.aggregate_script, cfg.get("aggregate_script")),
        "aggregate_sort_mem": pick(
            args.aggregate_sort_mem,
            cfg.get("aggregate_sort_mem", "8G"),
        ),
        "methscan_prepare_script": pick(
            args.methscan_prepare_script,
            cfg.get("methscan_prepare_script"),
        ),
        "methscan_pixi_manifest": pick(
            args.methscan_pixi_manifest,
            cfg.get("methscan_pixi_manifest"),
        ),
        "methscan_prepare_chunksize": pick(
            args.methscan_prepare_chunksize,
            cfg.get("methscan_prepare_chunksize"),
        ),
        "slurm_cfg_raw": slurm_cfg_raw,
        "slurm_partition": pick(args.slurm_partition, stage_slurm_cfg.get("partition")),
        "slurm_mem": pick(args.slurm_mem, stage_slurm_cfg.get("mem")),
        "slurm_cpus_per_task": pick(
            args.slurm_cpus_per_task, stage_slurm_cfg.get("cpus_per_task")
        ),
        "slurm_output": pick(args.slurm_output, stage_slurm_cfg.get("output")),
        "slurm_error": pick(args.slurm_error, stage_slurm_cfg.get("error")),
        "submit": args.submit,
        "dry_run": args.dry_run,
    }

    settings["work_root"] = settings["work_root"] or "work"
    settings["fastp_threads"] = (
        int(settings["fastp_threads"]) if settings["fastp_threads"] is not None else 8
    )
    if settings["fastp_threads"] <= 0:
        raise ValueError("fastp_threads must be > 0")

    if settings["number_of_split_parts"] is None:
        raise ValueError("number_of_split_parts is required")
    settings["number_of_split_parts"] = int(settings["number_of_split_parts"])
    if settings["number_of_split_parts"] <= 0:
        raise ValueError("number_of_split_parts must be > 0")

    if settings["biscuit_threads"] is not None:
        settings["biscuit_threads"] = int(settings["biscuit_threads"])
        if settings["biscuit_threads"] <= 0:
            raise ValueError("biscuit_threads must be > 0")
    else:
        settings["biscuit_threads"] = 32

    if settings["biscuit_batch_size"] is not None:
        settings["biscuit_batch_size"] = int(settings["biscuit_batch_size"])
        if settings["biscuit_batch_size"] <= 0:
            raise ValueError("biscuit_batch_size must be > 0")
    else:
        settings["biscuit_batch_size"] = 1

    settings["fastp_bin"] = normalize_executable_setting(
        settings["fastp_bin"], "fastp"
    )
    settings["biscuit_bin"] = normalize_executable_setting(
        settings["biscuit_bin"], "biscuit"
    )
    settings["sinto_bin"] = normalize_executable_setting(
        settings["sinto_bin"], "sinto"
    )
    settings["samtools_bin"] = normalize_executable_setting(
        settings["samtools_bin"], "samtools"
    )
    settings["call_bgzip_bin"] = normalize_executable_setting(
        settings["call_bgzip_bin"], "bgzip"
    )
    settings["call_tabix_bin"] = normalize_executable_setting(
        settings["call_tabix_bin"], "tabix"
    )
    settings["linker1"] = settings["linker1"]
    settings["linker2"] = settings["linker2"]
    settings["tn5"] = settings["tn5"]
    settings["linker_bc"] = (
        settings["linker_bc"]
        or settings["linker2"]
        or "ATTTATGTGTTTGAGAGGTTAGAGTATTTG"
    )
    settings["insert_left"] = (
        settings["insert_left"]
        or settings["tn5"]
        or "TATTGGTGTATGATTAGATGTGTATAAGAGATAG"
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

    settings["slurm_partition"] = settings["slurm_partition"] or "cpu"
    settings["slurm_mem"] = settings["slurm_mem"] or "16G"
    settings["slurm_cpus_per_task"] = settings["slurm_cpus_per_task"] or 8

    # Pool / split defaults.
    settings["samtools_threads"] = (
        int(settings["samtools_threads"])
        if settings["samtools_threads"] is not None
        else 4
    )
    if settings["samtools_threads"] <= 0:
        raise ValueError("samtools_threads must be > 0")
    settings["host_sort_mem"] = settings["host_sort_mem"] or "16G"

    settings["split_barcodes"] = settings["split_barcodes"] or cfg.get(
        "barcode1_whitelist"
    )
    settings["split_cb_tag"] = settings["split_cb_tag"] or "CB"
    settings["split_threads_read"] = (
        int(settings["split_threads_read"])
        if settings["split_threads_read"] is not None
        else 1
    )
    if settings["split_threads_read"] <= 0:
        raise ValueError("split_threads_read must be > 0")
    settings["split_threads_write"] = (
        int(settings["split_threads_write"])
        if settings["split_threads_write"] is not None
        else 1
    )
    if settings["split_threads_write"] <= 0:
        raise ValueError("split_threads_write must be > 0")
    settings["split_sort_jobs"] = (
        int(settings["split_sort_jobs"])
        if settings["split_sort_jobs"] is not None
        else 8
    )
    if settings["split_sort_jobs"] <= 0:
        raise ValueError("split_sort_jobs must be > 0")

    # Call defaults.
    settings["call_jobs"] = (
        int(settings["call_jobs"]) if settings["call_jobs"] is not None else 8
    )
    if settings["call_jobs"] <= 0:
        raise ValueError("call_jobs must be > 0")
    settings["call_mode"] = settings["call_mode"] or "all"
    if settings["call_mode"] not in {"all", "host", "spike"}:
        raise ValueError("call_mode must be one of: all, host, spike")

    settings["call_host_threads"] = (
        int(settings["call_host_threads"])
        if settings["call_host_threads"] is not None
        else int(settings["biscuit_threads"])
    )
    if settings["call_host_threads"] <= 0:
        raise ValueError("call_host_threads must be > 0")
    settings["call_spike_threads"] = (
        int(settings["call_spike_threads"])
        if settings["call_spike_threads"] is not None
        else int(settings["biscuit_threads"])
    )
    if settings["call_spike_threads"] <= 0:
        raise ValueError("call_spike_threads must be > 0")

    for trim_key in ("call_left_trimming", "call_right_trimming"):
        if settings[trim_key] is not None:
            settings[trim_key] = int(settings[trim_key])
            if settings[trim_key] < 0:
                raise ValueError(f"{trim_key} must be >= 0")

    settings["host_subsample_fraction"] = (
        float(settings["host_subsample_fraction"])
        if settings["host_subsample_fraction"] is not None
        else 0.1
    )
    if settings["host_subsample_fraction"] <= 0 or settings["host_subsample_fraction"] > 1:
        raise ValueError("host_subsample_fraction must be in (0, 1]")
    settings["host_subsample_seed"] = (
        int(settings["host_subsample_seed"])
        if settings["host_subsample_seed"] is not None
        else 11
    )
    if settings["host_subsample_seed"] < 0:
        raise ValueError("host_subsample_seed must be >= 0")
    settings["call_mito_chromosomes"] = settings["call_mito_chromosomes"] or "chrM"

    settings["saturation_script"] = settings["saturation_script"] or "scripts/saturation.py"
    settings["saturation_reads_threshold"] = (
        float(settings["saturation_reads_threshold"])
        if settings["saturation_reads_threshold"] is not None
        else 1_000_000.0
    )
    if settings["saturation_reads_threshold"] <= 0:
        raise ValueError("saturation_reads_threshold must be > 0")

    settings["summary_script"] = settings["summary_script"] or "scripts/summary.py"
    settings["aggregate_script"] = settings["aggregate_script"] or "scripts/aggregate.py"
    settings["methscan_prepare_script"] = (
        settings["methscan_prepare_script"] or "scripts/methscan_prepare.py"
    )
    if settings["methscan_prepare_chunksize"] is not None:
        settings["methscan_prepare_chunksize"] = int(settings["methscan_prepare_chunksize"])
    else:
        settings["methscan_prepare_chunksize"] = 10_000_000
    if settings["methscan_prepare_chunksize"] < 1:
        raise ValueError("methscan_prepare_chunksize must be >= 1")

    settings["mbias_script"] = settings["mbias_script"] or "scripts-emseq/mbias.py"
    settings["mbias_mode"] = settings["mbias_mode"] or "spike"
    if settings["mbias_mode"] not in {"all", "host", "spike"}:
        raise ValueError("mbias_mode must be one of: all, host, spike")
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

    # Step-specific slurm for split.
    settings["split_bams_slurm_partition"] = (
        pick(args.slurm_partition, split_bams_slurm_cfg.get("partition"))
        or settings["slurm_partition"]
    )
    settings["split_bams_slurm_mem"] = (
        pick(args.slurm_mem, split_bams_slurm_cfg.get("mem")) or settings["slurm_mem"]
    )
    settings["split_bams_slurm_cpus_per_task"] = (
        pick(args.slurm_cpus_per_task, split_bams_slurm_cfg.get("cpus_per_task"))
        or settings["slurm_cpus_per_task"]
    )

    settings["split_sort_slurm_partition"] = (
        pick(args.slurm_partition, split_sort_slurm_cfg.get("partition"))
        or settings["slurm_partition"]
    )
    settings["split_sort_slurm_mem"] = (
        pick(args.slurm_mem, split_sort_slurm_cfg.get("mem")) or settings["slurm_mem"]
    )
    settings["split_sort_slurm_cpus_per_task"] = (
        pick(args.slurm_cpus_per_task, split_sort_slurm_cfg.get("cpus_per_task"))
        or settings["slurm_cpus_per_task"]
    )

    # Step-specific slurm for call: slurm.call.host / slurm.call.spike.
    if isinstance(stage_slurm_cfg, dict) and ("host" in stage_slurm_cfg or "spike" in stage_slurm_cfg):
        call_host_cfg = stage_slurm_cfg.get("host") or {}
        call_spike_cfg = stage_slurm_cfg.get("spike") or {}
    else:
        call_host_cfg = stage_slurm_cfg
        call_spike_cfg = stage_slurm_cfg

    settings["call_host_slurm_partition"] = (
        pick(args.slurm_partition, call_host_cfg.get("partition"))
        or settings["slurm_partition"]
    )
    settings["call_host_slurm_mem"] = (
        pick(args.slurm_mem, call_host_cfg.get("mem")) or settings["slurm_mem"]
    )
    settings["call_host_slurm_cpus_per_task"] = (
        pick(args.slurm_cpus_per_task, call_host_cfg.get("cpus_per_task"))
        or settings["slurm_cpus_per_task"]
    )
    settings["call_spike_slurm_partition"] = (
        pick(args.slurm_partition, call_spike_cfg.get("partition"))
        or settings["slurm_partition"]
    )
    settings["call_spike_slurm_mem"] = (
        pick(args.slurm_mem, call_spike_cfg.get("mem")) or settings["slurm_mem"]
    )
    settings["call_spike_slurm_cpus_per_task"] = (
        pick(args.slurm_cpus_per_task, call_spike_cfg.get("cpus_per_task"))
        or settings["slurm_cpus_per_task"]
    )

    settings["mbias_host_slurm_partition"] = (
        pick(args.slurm_partition, mbias_host_slurm_cfg.get("partition"))
        or settings["slurm_partition"]
    )
    settings["mbias_host_slurm_mem"] = (
        pick(args.slurm_mem, mbias_host_slurm_cfg.get("mem")) or settings["slurm_mem"]
    )
    settings["mbias_host_slurm_cpus_per_task"] = (
        pick(args.slurm_cpus_per_task, mbias_host_slurm_cfg.get("cpus_per_task"))
        or settings["slurm_cpus_per_task"]
    )
    settings["mbias_spike_slurm_partition"] = (
        pick(args.slurm_partition, mbias_spike_slurm_cfg.get("partition"))
        or settings["slurm_partition"]
    )
    settings["mbias_spike_slurm_mem"] = (
        pick(args.slurm_mem, mbias_spike_slurm_cfg.get("mem")) or settings["slurm_mem"]
    )
    settings["mbias_spike_slurm_cpus_per_task"] = (
        pick(args.slurm_cpus_per_task, mbias_spike_slurm_cfg.get("cpus_per_task"))
        or settings["slurm_cpus_per_task"]
    )

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

    if stage == "all":
        for stage_name in STAGE_SEQUENCE:
            validate_required_for_stage(stage_name, settings)
    else:
        validate_required_for_stage(stage, settings)
    return settings


def normalize_spike_in_index(raw) -> list[str]:
    if raw is None:
        return []
    if isinstance(raw, list):
        return [str(item) for item in raw]
    if isinstance(raw, dict):
        items: list[str] = []
        for name, index in raw.items():
            items.append(f"{name}={index}")
        return items
    raise ValueError("spike_in_index must be either an object or an array")


def validate_mbias_settings(settings: dict) -> None:
    """Extra validation for mbias stage (not covered by STAGE_REQUIRED_FIELDS)."""
    mode = settings["mbias_mode"]
    if mode in ("all", "host") and not settings.get("call_reference_file"):
        raise ValueError(
            "mbias stage requires call_reference_file when mbias_mode is all or host"
        )
    if mode in ("all", "spike") and not parse_spike_names(settings["spike_in_index"]):
        raise ValueError(
            "mbias stage requires spike_in_index when mbias_mode is all or spike"
        )


def generate_local_script(command: str, output_path: Path) -> None:
    content = (
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n\n"
        f"{command}\n"
    )
    write_text(output_path, content)
    output_path.chmod(0o755)


def generate_slurm_script(
    command: str, output_path: Path, log_dir: Path, settings: dict
) -> None:
    log_dir.mkdir(parents=True, exist_ok=True)
    job_name = settings.get("job_name") or f"emseq_{settings['stage']}_{settings['sample_id']}"
    lines = [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --partition={settings['slurm_partition']}",
        f"#SBATCH --cpus-per-task={settings['slurm_cpus_per_task']}",
        f"#SBATCH --mem={settings['slurm_mem']}",
        f"#SBATCH --output={settings['slurm_output'].replace('%x', job_name)}",
        f"#SBATCH --error={settings['slurm_error'].replace('%x', job_name)}",
        "",
        "set -euo pipefail",
        "",
        command,
        "",
    ]
    content = "\n".join(lines)
    write_text(output_path, content)
    output_path.chmod(0o755)


def submit_script(path: Path, runner: str) -> None:
    if runner == "local":
        subprocess.run(["bash", str(path)], check=True)
    else:
        subprocess.run(["sbatch", str(path)], check=True)


def generate_split_slurm_submit_script(
    split_script_path: Path, sort_script_path: Path, output_path: Path
) -> None:
    """Create a helper that submits split_bams first, then sort via dependency."""
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


def parse_emseq_generated_paths(command_output: str) -> list[Path]:
    """Collect paths from child `make_cmd` stdout lines."""
    generated: list[Path] = []
    prefix = "[emseq.make_cmd] generated="
    for line in command_output.splitlines():
        if line.startswith(prefix):
            generated.append(Path(line[len(prefix) :].strip()))
    return generated


def build_stage_passthrough_args(argv: list[str]) -> list[str]:
    """Strip `--stage` / `--runner` / `--submit` / `--dry-run` for subprocess invocations."""
    passthrough: list[str] = []
    flags_without_value = {"--split-smoke", "--submit", "--dry-run"}
    index = 0
    while index < len(argv):
        token = argv[index]
        if token in {"--stage", "--runner"}:
            index += 2
            continue
        if token.startswith("--stage=") or token.startswith("--runner="):
            index += 1
            continue
        if token in {"--submit", "--dry-run"}:
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


def generate_local_driver_script(
    stage_scripts: list[tuple[str, list[Path]]], output_path: Path
) -> None:
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


def generate_emseq_slurm_driver_script(
    stage_scripts: list[tuple[str, list[Path]]],
    output_path: Path,
    log_dir: Path,
    driver_settings: dict,
) -> None:
    """Chain stage sbatch scripts with afterok; split uses internal two-step deps."""
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
        lines.append(f'echo "[emseq.run.sbatch] stage={stage_name}"')
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
    lines.append('echo "[emseq.run.sbatch] done final_dep=${prev_stage_deps}"')
    driver_command = "\n".join(lines)
    job_name = driver_settings["job_name"]
    slurm_settings = {
        "sample_id": driver_settings["sample_id"],
        "stage": "all",
        "slurm_partition": driver_settings["slurm_partition"],
        "slurm_mem": driver_settings["slurm_mem"],
        "slurm_cpus_per_task": driver_settings["slurm_cpus_per_task"],
        "slurm_output": driver_settings["slurm_output"],
        "slurm_error": driver_settings["slurm_error"],
        "job_name": job_name,
    }
    generate_slurm_script(driver_command, output_path, log_dir, slurm_settings)


def emseq_driver_slurm_settings(settings: dict) -> dict:
    """Use `slurm.summary` (or fallbacks) for the all-pipeline driver sbatch headers."""
    summary_cfg = select_stage_slurm_cfg(settings["slurm_cfg_raw"], "summary")
    sample_id = settings["sample_id"]
    work_root = Path(settings["work_root"])
    partition = summary_cfg.get("partition") or settings["slurm_partition"]
    mem = summary_cfg.get("mem") or settings["slurm_mem"]
    cpus = summary_cfg.get("cpus_per_task") or settings["slurm_cpus_per_task"]
    job_name = f"emseq_all_driver_{sample_id}"
    out_tpl = summary_cfg.get("output") or str(
        work_root / sample_id / "logs" / "all_driver_%x_%j.out"
    )
    err_tpl = summary_cfg.get("error") or str(
        work_root / sample_id / "logs" / "all_driver_%x_%j.err"
    )
    return {
        "sample_id": sample_id,
        "job_name": job_name,
        "slurm_partition": partition,
        "slurm_mem": mem,
        "slurm_cpus_per_task": cpus,
        "slurm_output": out_tpl.replace("%x", job_name),
        "slurm_error": err_tpl.replace("%x", job_name),
    }


def main() -> int:
    args = parse_args()
    settings = resolve_settings(args)
    settings["spike_in_index"] = normalize_spike_in_index(settings["spike_in_index"])
    if settings["stage"] in ("mbias", "all"):
        validate_mbias_settings(settings)

    sample_work = Path(settings["work_root"]) / settings["sample_id"]
    command_dir = sample_work / "commands"
    log_dir = sample_work / "logs"

    stage = settings["stage"]

    if stage == "all":
        passthrough_args = build_stage_passthrough_args(sys.argv[1:])
        make_cmd_path = Path(__file__).resolve()
        stage_scripts: list[tuple[str, list[Path]]] = []
        for stage_name in STAGE_SEQUENCE:
            stage_argv = [
                sys.executable,
                str(make_cmd_path),
                *passthrough_args,
                "--runner",
                settings["runner"],
                "--stage",
                stage_name,
            ]
            if settings["dry_run"]:
                stage_argv.append("--dry-run")
            completed = subprocess.run(
                stage_argv, check=True, capture_output=True, text=True
            )
            if completed.stdout:
                print(completed.stdout, end="")
            if completed.stderr:
                print(completed.stderr, end="", file=sys.stderr)
            combined_out = (completed.stdout or "") + (completed.stderr or "")
            stage_scripts.append(
                (stage_name, parse_emseq_generated_paths(combined_out))
            )

        driver_path: Path
        if settings["runner"] == "local":
            driver_path = command_dir / "run.sh"
            print(f"[emseq.make_cmd] script={driver_path}")
            if not settings["dry_run"]:
                generate_local_driver_script(stage_scripts, driver_path)
        else:
            driver_path = command_dir / "run.sbatch"
            print(f"[emseq.make_cmd] script={driver_path}")
            if not settings["dry_run"]:
                driver_slurm = emseq_driver_slurm_settings(settings)
                generate_emseq_slurm_driver_script(
                    stage_scripts, driver_path, log_dir, driver_slurm
                )

        if not settings["dry_run"] and driver_path.exists():
            print(f"[emseq.make_cmd] generated={driver_path}")
        if settings["submit"] and not settings["dry_run"]:
            if settings["runner"] == "slurm":
                subprocess.run(["bash", str(driver_path)], check=True)
                print("[emseq.make_cmd] submitted_driver=1")
                print("[emseq.make_cmd] submit_mode=client_side_sbatch_dag")
            else:
                submit_script(driver_path, settings["runner"])
                print("[emseq.make_cmd] submitted_driver=1")

        print("[emseq.make_cmd] stage=all helper generation complete")
        return 0

    if stage == "fastp_split":
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

        print(f"[emseq.make_cmd] runner={settings['runner']}")
        print(f"[emseq.make_cmd] stage={stage}")
        print(f"[emseq.make_cmd] sample_id={settings['sample_id']}")
        print(f"[emseq.make_cmd] script={script_path}")
        print(f"[emseq.make_cmd] command={command}")

        if settings["dry_run"]:
            return 0

        if settings["runner"] == "local":
            generate_local_script(command, script_path)
        else:
            slurm_settings = dict(settings)
            slurm_settings["stage"] = stage
            generate_slurm_script(command, script_path, log_dir, slurm_settings)

        print(f"[emseq.make_cmd] generated={script_path}")

        if settings["submit"]:
            submit_script(script_path, settings["runner"])
            print("[emseq.make_cmd] submitted_count=1")
    elif stage == "demux_extract_bc":
        command_args = argparse.Namespace(
            barcode1_whitelist=settings["barcode1_whitelist"],
            barcode2_whitelist=settings["barcode2_whitelist"],
            linker_bc=settings["linker_bc"],
            insert_left=settings["insert_left"],
            linker_edit_distance=settings["linker_edit_distance"],
            barcode_hamming_distance=settings["barcode_hamming_distance"],
            gzip_level=settings["gzip_level"],
        )

        print(f"[emseq.make_cmd] runner={settings['runner']}")
        print(f"[emseq.make_cmd] stage={stage}")
        print(f"[emseq.make_cmd] sample_id={settings['sample_id']}")

        if settings["runner"] == "local":
            script_path = command_dir / "02_demux_extract_bc.sh"
            command = build_demux_local_batch_command(command_args, sample_work)
            print(f"[emseq.make_cmd] script={script_path}")
            print(f"[emseq.make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            print(f"[emseq.make_cmd] generated={script_path}")
            if settings["submit"]:
                submit_script(script_path, settings["runner"])
                print("[emseq.make_cmd] submitted_count=1")
        else:
            chunks = build_demux_chunks_from_config(
                sample_work, settings["number_of_split_parts"]
            )
            print(f"[emseq.make_cmd] chunk_count={len(chunks)}")
            generated_scripts: list[Path] = []
            for chunk, r1_path, r2_path, out_prefix in chunks:
                base_name = f"02_demux_extract_bc_{chunk}"
                script_path = command_dir / f"{base_name}.sbatch"
                command = build_demux_chunk_command(
                    command_args, r1_path, r2_path, out_prefix
                )
                job_name = f"emseq_demux_{settings['sample_id']}_{chunk}"
                chunk_output = (settings["slurm_output"] or "").replace(
                    "%x", job_name
                )
                chunk_error = (settings["slurm_error"] or "").replace(
                    "%x", job_name
                )
                slurm_settings = {
                    "sample_id": settings["sample_id"],
                    "stage": stage,
                    "slurm_partition": settings["slurm_partition"],
                    "slurm_mem": settings["slurm_mem"],
                    "slurm_cpus_per_task": settings["slurm_cpus_per_task"],
                    "slurm_output": chunk_output,
                    "slurm_error": chunk_error,
                    "job_name": job_name,
                }
                print(f"[emseq.make_cmd] script={script_path}")
                print(f"[emseq.make_cmd] command={command}")
                if not settings["dry_run"]:
                    generate_slurm_script(command, script_path, log_dir, slurm_settings)
                    print(f"[emseq.make_cmd] generated={script_path}")
                generated_scripts.append(script_path)

            if settings["submit"] and not settings["dry_run"]:
                for script_path in generated_scripts:
                    submit_script(script_path, settings["runner"])
                print(f"[emseq.make_cmd] submitted_count={len(generated_scripts)}")
    elif stage == "align":
        command_args = argparse.Namespace(
            biscuit_reference=settings["biscuit_reference"],
            biscuit_threads=settings["biscuit_threads"],
            biscuit_batch_size=settings["biscuit_batch_size"],
            biscuit_bin=settings["biscuit_bin"],
            sinto_bin=settings["sinto_bin"],
            samtools_bin=settings["samtools_bin"],
            spike_in_index=settings["spike_in_index"],
        )

        print(f"[emseq.make_cmd] runner={settings['runner']}")
        print(f"[emseq.make_cmd] stage={stage}")
        print(f"[emseq.make_cmd] sample_id={settings['sample_id']}")

        if settings["runner"] == "local":
            base_name = "03_align"
            script_path = command_dir / f"{base_name}.sh"
            script = quoted(
                [
                    sys.executable,
                    "scripts-emseq/aligner.py",
                    "--work-path",
                    str(sample_work),
                    "--biscuit-reference",
                    command_args.biscuit_reference,
                    "--biscuit-threads",
                    str(command_args.biscuit_threads),
                    "--biscuit-batch-size",
                    str(command_args.biscuit_batch_size),
                    "--biscuit-bin",
                    command_args.biscuit_bin,
                    "--sinto-bin",
                    command_args.sinto_bin,
                    "--samtools-bin",
                    command_args.samtools_bin,
                ]
                + [
                    item
                    for value in command_args.spike_in_index
                    for item in ["--spike-in-index", value]
                ]
            )
            print(f"[emseq.make_cmd] script={script_path}")
            print(f"[emseq.make_cmd] command={script}")
            if settings["dry_run"]:
                return 0
            generate_local_script(script, script_path)
            print(f"[emseq.make_cmd] generated={script_path}")
            if settings["submit"]:
                submit_script(script_path, settings["runner"])
                print("[emseq.make_cmd] submitted_count=1")
        else:
            chunks = build_chunk_names(settings["number_of_split_parts"])
            print(f"[emseq.make_cmd] chunk_count={len(chunks)}")
            generated_scripts: list[Path] = []
            for chunk in chunks:
                base_name = f"03_align_{chunk}"
                script_path = command_dir / f"{base_name}.sbatch"
                command = quoted(
                    [
                        sys.executable,
                        "scripts-emseq/aligner.py",
                        "--work-path",
                        str(sample_work),
                        "--chunk",
                        chunk,
                        "--biscuit-reference",
                        command_args.biscuit_reference,
                        "--biscuit-threads",
                        str(command_args.biscuit_threads),
                        "--biscuit-batch-size",
                        str(command_args.biscuit_batch_size),
                        "--biscuit-bin",
                        command_args.biscuit_bin,
                        "--sinto-bin",
                        command_args.sinto_bin,
                        "--samtools-bin",
                        command_args.samtools_bin,
                    ]
                    + [
                        item
                        for value in command_args.spike_in_index
                        for item in ["--spike-in-index", value]
                    ]
                )
                job_name = f"emseq_align_{settings['sample_id']}_{chunk}"
                chunk_output = (settings["slurm_output"] or "").replace(
                    "%x", job_name
                )
                chunk_error = (settings["slurm_error"] or "").replace("%x", job_name)
                slurm_settings = {
                    "sample_id": settings["sample_id"],
                    "stage": stage,
                    "slurm_partition": settings["slurm_partition"],
                    "slurm_mem": settings["slurm_mem"],
                    "slurm_cpus_per_task": settings["slurm_cpus_per_task"],
                    "slurm_output": chunk_output,
                    "slurm_error": chunk_error,
                    "job_name": job_name,
                }
                print(f"[emseq.make_cmd] script={script_path}")
                print(f"[emseq.make_cmd] command={command}")
                if not settings["dry_run"]:
                    generate_slurm_script(command, script_path, log_dir, slurm_settings)
                    print(f"[emseq.make_cmd] generated={script_path}")
                generated_scripts.append(script_path)

            if settings["submit"] and not settings["dry_run"]:
                for script_path in generated_scripts:
                    submit_script(script_path, settings["runner"])
                print(f"[emseq.make_cmd] submitted_count={len(generated_scripts)}")
    elif stage == "pool":
        spike_names = parse_spike_names(settings["spike_in_index"])
        mode = "all" if spike_names else "host"
        command_args = argparse.Namespace(
            samtools_bin=settings["samtools_bin"],
            samtools_threads=settings["samtools_threads"],
            host_sort_mem=settings["host_sort_mem"],
            spike_in_index=settings["spike_in_index"],
        )

        print(f"[emseq.make_cmd] runner={settings['runner']}")
        print(f"[emseq.make_cmd] stage={stage}")
        print(f"[emseq.make_cmd] sample_id={settings['sample_id']}")
        print(f"[emseq.make_cmd] spike_in_count={len(spike_names)}")

        if settings["runner"] == "local":
            script_path = command_dir / "04_pool.sh"
            command = build_pool_command(command_args, sample_work, mode)
            print(f"[emseq.make_cmd] script={script_path}")
            print(f"[emseq.make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            print(f"[emseq.make_cmd] generated={script_path}")
            if settings["submit"]:
                submit_script(script_path, settings["runner"])
                print("[emseq.make_cmd] submitted_count=1")
        else:
            generated_scripts: list[Path] = []

            if spike_names:
                spike_script_path = command_dir / "04_pool_spike.sbatch"
                spike_job_name = f"emseq_pool_spike_{settings['sample_id']}"
                spike_command = build_pool_command(command_args, sample_work, "spike")
                spike_output = (settings["slurm_output"] or "").replace(
                    "%x", spike_job_name
                )
                spike_error = (settings["slurm_error"] or "").replace(
                    "%x", spike_job_name
                )
                spike_slurm_settings = {
                    "sample_id": settings["sample_id"],
                    "stage": stage,
                    "slurm_partition": settings["slurm_partition"],
                    "slurm_mem": settings["slurm_mem"],
                    "slurm_cpus_per_task": settings["slurm_cpus_per_task"],
                    "slurm_output": spike_output,
                    "slurm_error": spike_error,
                    "job_name": spike_job_name,
                }
                print(f"[emseq.make_cmd] script={spike_script_path}")
                print(f"[emseq.make_cmd] command={spike_command}")
                if not settings["dry_run"]:
                    generate_slurm_script(
                        spike_command, spike_script_path, log_dir, spike_slurm_settings
                    )
                    print(f"[emseq.make_cmd] generated={spike_script_path}")
                generated_scripts.append(spike_script_path)

            host_script_path = command_dir / "04_pool_host.sbatch"
            host_job_name = f"emseq_pool_host_{settings['sample_id']}"
            host_command = build_pool_command(command_args, sample_work, "host")
            host_output = (settings["slurm_output"] or "").replace(
                "%x", host_job_name
            )
            host_error = (settings["slurm_error"] or "").replace(
                "%x", host_job_name
            )
            host_slurm_settings = {
                "sample_id": settings["sample_id"],
                "stage": stage,
                "slurm_partition": settings["slurm_partition"],
                "slurm_mem": settings["slurm_mem"],
                "slurm_cpus_per_task": settings["slurm_cpus_per_task"],
                "slurm_output": host_output,
                "slurm_error": host_error,
                "job_name": host_job_name,
            }
            print(f"[emseq.make_cmd] script={host_script_path}")
            print(f"[emseq.make_cmd] command={host_command}")
            if not settings["dry_run"]:
                generate_slurm_script(
                    host_command, host_script_path, log_dir, host_slurm_settings
                )
                print(f"[emseq.make_cmd] generated={host_script_path}")
            generated_scripts.append(host_script_path)

            if settings["submit"] and not settings["dry_run"]:
                for script_path in generated_scripts:
                    submit_script(script_path, settings["runner"])
                print(
                    f"[emseq.make_cmd] submitted_count={len(generated_scripts)}"
                )
    elif stage == "split":
        command_args = argparse.Namespace(
            split_barcodes=settings["split_barcodes"],
            split_cb_tag=settings["split_cb_tag"],
            split_threads_read=settings["split_threads_read"],
            split_threads_write=settings["split_threads_write"],
            split_smoke=settings["split_smoke"],
            samtools_bin=settings["samtools_bin"],
            split_sort_jobs=settings["split_sort_jobs"],
        )

        print(f"[emseq.make_cmd] runner={settings['runner']}")
        print(f"[emseq.make_cmd] stage={stage}")
        print(f"[emseq.make_cmd] sample_id={settings['sample_id']}")

        if settings["runner"] == "local":
            script_path = command_dir / "05_split.sh"
            split_command = build_split_command(command_args, sample_work)
            sort_command = build_split_sort_command(command_args, sample_work)
            command = f"{split_command}\n{sort_command}"
            print(f"[emseq.make_cmd] script={script_path}")
            print(f"[emseq.make_cmd] command={split_command}")
            print(f"[emseq.make_cmd] command={sort_command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            print(f"[emseq.make_cmd] generated={script_path}")
            if settings["submit"]:
                submit_script(script_path, settings["runner"])
                print("[emseq.make_cmd] submitted_count=1")
        else:
            split_script_path = command_dir / "05_split_bams.sbatch"
            sort_script_path = command_dir / "05_split_sort.sbatch"
            split_submit_helper_path = command_dir / "05_split_submit.sh"
            split_command = build_split_command(command_args, sample_work)
            sort_command = build_split_sort_command(command_args, sample_work)

            split_job_name = f"emseq_split_bams_{settings['sample_id']}"
            sort_job_name = f"emseq_split_sort_{settings['sample_id']}"

            split_output = (settings["slurm_output"] or "").replace(
                "%x", split_job_name
            )
            split_error = (settings["slurm_error"] or "").replace(
                "%x", split_job_name
            )
            sort_output = (settings["slurm_output"] or "").replace(
                "%x", sort_job_name
            )
            sort_error = (settings["slurm_error"] or "").replace(
                "%x", sort_job_name
            )

            split_slurm_settings = {
                "sample_id": settings["sample_id"],
                "stage": stage,
                "slurm_partition": settings["split_bams_slurm_partition"],
                "slurm_mem": settings["split_bams_slurm_mem"],
                "slurm_cpus_per_task": settings["split_bams_slurm_cpus_per_task"],
                "slurm_output": split_output,
                "slurm_error": split_error,
                "job_name": split_job_name,
            }
            sort_slurm_settings = {
                "sample_id": settings["sample_id"],
                "stage": stage,
                "slurm_partition": settings["split_sort_slurm_partition"],
                "slurm_mem": settings["split_sort_slurm_mem"],
                "slurm_cpus_per_task": settings["split_sort_slurm_cpus_per_task"],
                "slurm_output": sort_output,
                "slurm_error": sort_error,
                "job_name": sort_job_name,
            }

            print(f"[emseq.make_cmd] script={split_script_path}")
            print(f"[emseq.make_cmd] command={split_command}")
            print(f"[emseq.make_cmd] script={sort_script_path}")
            print(f"[emseq.make_cmd] command={sort_command}")
            print(f"[emseq.make_cmd] helper={split_submit_helper_path}")

            if not settings["dry_run"]:
                generate_slurm_script(
                    split_command, split_script_path, log_dir, split_slurm_settings
                )
                print(f"[emseq.make_cmd] generated={split_script_path}")
                generate_slurm_script(
                    sort_command, sort_script_path, log_dir, sort_slurm_settings
                )
                print(f"[emseq.make_cmd] generated={sort_script_path}")
                generate_split_slurm_submit_script(
                    split_script_path,
                    sort_script_path,
                    split_submit_helper_path,
                )

                if settings["submit"]:
                    subprocess.run(
                        ["bash", str(split_submit_helper_path)], check=True
                    )
                    print("[emseq.make_cmd] submitted_count=1")
    elif stage == "mbias":
        command_args = argparse.Namespace(
            samtools_bin=settings["samtools_bin"],
            samtools_threads=settings["samtools_threads"],
            mbias_host_subsample_fraction=settings["mbias_host_subsample_fraction"],
            mbias_max_cycle=settings["mbias_max_cycle"],
            mbias_min_mapping_quality=settings["mbias_min_mapping_quality"],
            mbias_script=settings["mbias_script"],
            call_reference_file=settings["call_reference_file"],
            spike_in_index=settings["spike_in_index"],
        )

        print(f"[emseq.make_cmd] runner={settings['runner']}")
        print(f"[emseq.make_cmd] stage={stage}")
        print(f"[emseq.make_cmd] sample_id={settings['sample_id']}")
        print(f"[emseq.make_cmd] mbias_mode={settings['mbias_mode']}")

        if settings["runner"] == "local":
            script_path = command_dir / "06_mbias.sh"
            command = build_mbias_command(
                command_args,
                sample_work,
                settings["mbias_mode"],
            )
            print(f"[emseq.make_cmd] script={script_path}")
            print(f"[emseq.make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            print(f"[emseq.make_cmd] generated={script_path}")
            if settings["submit"]:
                submit_script(script_path, settings["runner"])
                print("[emseq.make_cmd] submitted_count=1")
        else:
            mbias_mode = settings["mbias_mode"]
            include_host = mbias_mode in ("all", "host")
            include_spike = mbias_mode in ("all", "spike")
            generated_scripts: list[Path] = []

            if include_host:
                host_script_path = command_dir / "06_mbias_host.sbatch"
                host_job_name = f"emseq_mbias_host_{settings['sample_id']}"
                host_command = build_mbias_command(
                    command_args, sample_work, "host"
                )
                host_output = (settings["slurm_output"] or "").replace(
                    "%x", host_job_name
                )
                host_error = (settings["slurm_error"] or "").replace(
                    "%x", host_job_name
                )
                host_slurm_settings = {
                    "sample_id": settings["sample_id"],
                    "stage": stage,
                    "slurm_partition": settings["mbias_host_slurm_partition"],
                    "slurm_mem": settings["mbias_host_slurm_mem"],
                    "slurm_cpus_per_task": settings["mbias_host_slurm_cpus_per_task"],
                    "slurm_output": host_output,
                    "slurm_error": host_error,
                    "job_name": host_job_name,
                }
                print(f"[emseq.make_cmd] script={host_script_path}")
                print(f"[emseq.make_cmd] command={host_command}")
                if not settings["dry_run"]:
                    generate_slurm_script(
                        host_command, host_script_path, log_dir, host_slurm_settings
                    )
                    print(f"[emseq.make_cmd] generated={host_script_path}")
                generated_scripts.append(host_script_path)

            if include_spike:
                for spike_name in parse_spike_names(settings["spike_in_index"]):
                    spike_script_path = (
                        command_dir / f"06_mbias_spike_{spike_name}.sbatch"
                    )
                    spike_job_name = (
                        f"emseq_mbias_spike_{settings['sample_id']}_{spike_name}"
                    )
                    spike_command = build_mbias_command(
                        command_args,
                        sample_work,
                        "spike",
                        spike_name=spike_name,
                    )
                    spike_output = (settings["slurm_output"] or "").replace(
                        "%x", spike_job_name
                    )
                    spike_error = (settings["slurm_error"] or "").replace(
                        "%x", spike_job_name
                    )
                    spike_slurm_settings = {
                        "sample_id": settings["sample_id"],
                        "stage": stage,
                        "slurm_partition": settings["mbias_spike_slurm_partition"],
                        "slurm_mem": settings["mbias_spike_slurm_mem"],
                        "slurm_cpus_per_task": settings["mbias_spike_slurm_cpus_per_task"],
                        "slurm_output": spike_output,
                        "slurm_error": spike_error,
                        "job_name": spike_job_name,
                    }
                    print(f"[emseq.make_cmd] script={spike_script_path}")
                    print(f"[emseq.make_cmd] command={spike_command}")
                    if not settings["dry_run"]:
                        generate_slurm_script(
                            spike_command,
                            spike_script_path,
                            log_dir,
                            spike_slurm_settings,
                        )
                        print(f"[emseq.make_cmd] generated={spike_script_path}")
                    generated_scripts.append(spike_script_path)

            if settings["submit"] and not settings["dry_run"]:
                for script_path in generated_scripts:
                    submit_script(script_path, settings["runner"])
                print(
                    f"[emseq.make_cmd] submitted_count={len(generated_scripts)}"
                )
    elif stage == "call":
        spike_names = parse_spike_names(settings["spike_in_index"])
        effective_mode = settings["call_mode"]
        if effective_mode in ("all", "spike") and not spike_names:
            print("[emseq.make_cmd] call spike_in_index empty: run host only")
            effective_mode = "host"

        include_host = effective_mode in ("all", "host")
        include_spike = effective_mode in ("all", "spike")

        command_args = argparse.Namespace(
            call_reference_file=settings["call_reference_file"],
            call_jobs=settings["call_jobs"],
            call_host_threads=settings["call_host_threads"],
            call_spike_threads=settings["call_spike_threads"],
            call_bgzip_bin=settings["call_bgzip_bin"],
            call_tabix_bin=settings["call_tabix_bin"],
            call_mito_chromosomes=settings["call_mito_chromosomes"],
            call_left_trimming=settings["call_left_trimming"],
            call_right_trimming=settings["call_right_trimming"],
            biscuit_bin=settings["biscuit_bin"],
            spike_in_index=settings["spike_in_index"],
        )

        print(f"[emseq.make_cmd] runner={settings['runner']}")
        print(f"[emseq.make_cmd] stage={stage}")
        print(f"[emseq.make_cmd] sample_id={settings['sample_id']}")
        print(f"[emseq.make_cmd] call_mode={effective_mode}")
        print(f"[emseq.make_cmd] spike_in_count={len(spike_names)}")

        if settings["runner"] == "local":
            script_path = command_dir / "08_call.sh"
            spike_references = (
                settings["spike_in_index"] if include_spike else []
            )
            command = build_call_command(
                command_args,
                sample_work,
                effective_mode,
                spike_references=spike_references,
            )
            print(f"[emseq.make_cmd] script={script_path}")
            print(f"[emseq.make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            print(f"[emseq.make_cmd] generated={script_path}")
            if settings["submit"]:
                submit_script(script_path, settings["runner"])
                print("[emseq.make_cmd] submitted_count=1")
        else:
            generated_scripts: list[Path] = []

            if include_host:
                host_script_path = command_dir / "08_call_host.sbatch"
                host_job_name = f"emseq_call_host_{settings['sample_id']}"
                host_output = (settings["slurm_output"] or "").replace(
                    "%x", host_job_name
                )
                host_error = (settings["slurm_error"] or "").replace(
                    "%x", host_job_name
                )
                host_command = build_call_command(
                    command_args,
                    sample_work,
                    "host",
                    spike_references=[],
                )
                host_slurm_settings = {
                    "sample_id": settings["sample_id"],
                    "stage": stage,
                    "slurm_partition": settings["call_host_slurm_partition"],
                    "slurm_mem": settings["call_host_slurm_mem"],
                    "slurm_cpus_per_task": settings["call_host_slurm_cpus_per_task"],
                    "slurm_output": host_output,
                    "slurm_error": host_error,
                    "job_name": host_job_name,
                }
                print(f"[emseq.make_cmd] script={host_script_path}")
                print(f"[emseq.make_cmd] command={host_command}")
                if not settings["dry_run"]:
                    generate_slurm_script(
                        host_command, host_script_path, log_dir, host_slurm_settings
                    )
                    print(f"[emseq.make_cmd] generated={host_script_path}")
                generated_scripts.append(host_script_path)

            if include_spike:
                for spike_name in spike_names:
                    spike_script_path = (
                        command_dir / f"08_call_spike_{spike_name}.sbatch"
                    )
                    spike_job_name = (
                        f"emseq_call_spike_{settings['sample_id']}_{spike_name}"
                    )
                    spike_output = (settings["slurm_output"] or "").replace(
                        "%x", spike_job_name
                    )
                    spike_error = (settings["slurm_error"] or "").replace(
                        "%x", spike_job_name
                    )
                    spike_reference_items: list[str] = [
                        item
                        for item in settings["spike_in_index"]
                        if item.split("=", 1)[0] == spike_name
                    ]
                    spike_command = build_call_command(
                        command_args,
                        sample_work,
                        "spike",
                        spike_name=spike_name,
                        spike_references=spike_reference_items,
                    )
                    spike_slurm_settings = {
                        "sample_id": settings["sample_id"],
                        "stage": stage,
                        "slurm_partition": settings["call_spike_slurm_partition"],
                        "slurm_mem": settings["call_spike_slurm_mem"],
                        "slurm_cpus_per_task": settings["call_spike_slurm_cpus_per_task"],
                        "slurm_output": spike_output,
                        "slurm_error": spike_error,
                        "job_name": spike_job_name,
                    }
                    print(f"[emseq.make_cmd] script={spike_script_path}")
                    print(f"[emseq.make_cmd] command={spike_command}")
                    if not settings["dry_run"]:
                        generate_slurm_script(
                            spike_command,
                            spike_script_path,
                            log_dir,
                            spike_slurm_settings,
                        )
                        print(f"[emseq.make_cmd] generated={spike_script_path}")
                    generated_scripts.append(spike_script_path)

            if settings["submit"] and not settings["dry_run"]:
                for script_path in generated_scripts:
                    submit_script(script_path, settings["runner"])
                print(
                    f"[emseq.make_cmd] submitted_count={len(generated_scripts)}"
                )
    elif stage == "saturation":
        command_args = argparse.Namespace(
            saturation_script=settings["saturation_script"],
            saturation_reads_threshold=settings["saturation_reads_threshold"],
        )
        command = build_saturation_command(command_args, sample_work)
        if settings["runner"] == "local":
            script_path = command_dir / "09_saturation.sh"
            print(f"[emseq.make_cmd] runner={settings['runner']}")
            print(f"[emseq.make_cmd] stage={stage}")
            print(f"[emseq.make_cmd] sample_id={settings['sample_id']}")
            print(f"[emseq.make_cmd] script={script_path}")
            print(f"[emseq.make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            print(f"[emseq.make_cmd] generated={script_path}")
            if settings["submit"]:
                submit_script(script_path, settings["runner"])
                print("[emseq.make_cmd] submitted_count=1")
        else:
            script_path = command_dir / "09_saturation.sbatch"
            job_name = f"emseq_saturation_{settings['sample_id']}"
            sat_output = (settings["slurm_output"] or "").replace("%x", job_name)
            sat_error = (settings["slurm_error"] or "").replace("%x", job_name)
            slurm_settings = {
                "sample_id": settings["sample_id"],
                "stage": stage,
                "slurm_partition": settings["slurm_partition"],
                "slurm_mem": settings["slurm_mem"],
                "slurm_cpus_per_task": settings["slurm_cpus_per_task"],
                "slurm_output": sat_output,
                "slurm_error": sat_error,
                "job_name": job_name,
            }
            print(f"[emseq.make_cmd] runner={settings['runner']}")
            print(f"[emseq.make_cmd] stage={stage}")
            print(f"[emseq.make_cmd] sample_id={settings['sample_id']}")
            print(f"[emseq.make_cmd] script={script_path}")
            print(f"[emseq.make_cmd] command={command}")
            if not settings["dry_run"]:
                generate_slurm_script(command, script_path, log_dir, slurm_settings)
                print(f"[emseq.make_cmd] generated={script_path}")
            if settings["submit"] and not settings["dry_run"]:
                submit_script(script_path, settings["runner"])
                print("[emseq.make_cmd] submitted_count=1")
    elif stage == "summary":
        command_args = argparse.Namespace(
            summary_script=settings["summary_script"],
            spike_in_index=settings["spike_in_index"],
        )
        command = build_summary_command(command_args, sample_work)
        if settings["runner"] == "local":
            script_path = command_dir / "10_summary.sh"
            print(f"[emseq.make_cmd] runner={settings['runner']}")
            print(f"[emseq.make_cmd] stage={stage}")
            print(f"[emseq.make_cmd] sample_id={settings['sample_id']}")
            print(f"[emseq.make_cmd] script={script_path}")
            print(f"[emseq.make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            print(f"[emseq.make_cmd] generated={script_path}")
            if settings["submit"]:
                submit_script(script_path, settings["runner"])
                print("[emseq.make_cmd] submitted_count=1")
        else:
            script_path = command_dir / "10_summary.sbatch"
            job_name = f"emseq_summary_{settings['sample_id']}"
            summary_output = (settings["slurm_output"] or "").replace("%x", job_name)
            summary_error = (settings["slurm_error"] or "").replace("%x", job_name)
            slurm_settings = {
                "sample_id": settings["sample_id"],
                "stage": stage,
                "slurm_partition": settings["slurm_partition"],
                "slurm_mem": settings["slurm_mem"],
                "slurm_cpus_per_task": settings["slurm_cpus_per_task"],
                "slurm_output": summary_output,
                "slurm_error": summary_error,
                "job_name": job_name,
            }
            print(f"[emseq.make_cmd] runner={settings['runner']}")
            print(f"[emseq.make_cmd] stage={stage}")
            print(f"[emseq.make_cmd] sample_id={settings['sample_id']}")
            print(f"[emseq.make_cmd] script={script_path}")
            print(f"[emseq.make_cmd] command={command}")
            if not settings["dry_run"]:
                generate_slurm_script(command, script_path, log_dir, slurm_settings)
                print(f"[emseq.make_cmd] generated={script_path}")
            if settings["submit"] and not settings["dry_run"]:
                submit_script(script_path, settings["runner"])
                print("[emseq.make_cmd] submitted_count=1")
    elif stage == "aggregate":
        command_args = argparse.Namespace(
            aggregate_script=settings["aggregate_script"],
            sort_mem=settings["aggregate_sort_mem"],
        )
        command = build_aggregate_command(command_args, sample_work)
        if settings["runner"] == "local":
            script_path = command_dir / "11_aggregate.sh"
            print(f"[emseq.make_cmd] runner={settings['runner']}")
            print(f"[emseq.make_cmd] stage={stage}")
            print(f"[emseq.make_cmd] sample_id={settings['sample_id']}")
            print(f"[emseq.make_cmd] script={script_path}")
            print(f"[emseq.make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            print(f"[emseq.make_cmd] generated={script_path}")
            if settings["submit"]:
                submit_script(script_path, settings["runner"])
                print("[emseq.make_cmd] submitted_count=1")
        else:
            script_path = command_dir / "11_aggregate.sbatch"
            job_name = f"emseq_aggregate_{settings['sample_id']}"
            aggregate_output = (settings["slurm_output"] or "").replace("%x", job_name)
            aggregate_error = (settings["slurm_error"] or "").replace("%x", job_name)
            slurm_settings = {
                "sample_id": settings["sample_id"],
                "stage": stage,
                "slurm_partition": settings["slurm_partition"],
                "slurm_mem": settings["slurm_mem"],
                "slurm_cpus_per_task": settings["slurm_cpus_per_task"],
                "slurm_output": aggregate_output,
                "slurm_error": aggregate_error,
                "job_name": job_name,
            }
            print(f"[emseq.make_cmd] runner={settings['runner']}")
            print(f"[emseq.make_cmd] stage={stage}")
            print(f"[emseq.make_cmd] sample_id={settings['sample_id']}")
            print(f"[emseq.make_cmd] script={script_path}")
            print(f"[emseq.make_cmd] command={command}")
            if not settings["dry_run"]:
                generate_slurm_script(command, script_path, log_dir, slurm_settings)
                print(f"[emseq.make_cmd] generated={script_path}")
            if settings["submit"] and not settings["dry_run"]:
                submit_script(script_path, settings["runner"])
                print("[emseq.make_cmd] submitted_count=1")
    elif stage == "methscan_prepare":
        command_args = argparse.Namespace(
            methscan_prepare_script=settings["methscan_prepare_script"],
            methscan_pixi_manifest=settings.get("methscan_pixi_manifest"),
            methscan_prepare_chunksize=settings["methscan_prepare_chunksize"],
        )
        command = build_methscan_prepare_command(command_args, sample_work)
        if settings["runner"] == "local":
            script_path = command_dir / "12_methscan_prepare.sh"
            print(f"[emseq.make_cmd] runner={settings['runner']}")
            print(f"[emseq.make_cmd] stage={stage}")
            print(f"[emseq.make_cmd] sample_id={settings['sample_id']}")
            print(f"[emseq.make_cmd] script={script_path}")
            print(f"[emseq.make_cmd] command={command}")
            if settings["dry_run"]:
                return 0
            generate_local_script(command, script_path)
            print(f"[emseq.make_cmd] generated={script_path}")
            if settings["submit"]:
                submit_script(script_path, settings["runner"])
                print("[emseq.make_cmd] submitted_count=1")
        else:
            script_path = command_dir / "12_methscan_prepare.sbatch"
            job_name = f"emseq_methscan_prepare_{settings['sample_id']}"
            methscan_output = (settings["slurm_output"] or "").replace("%x", job_name)
            methscan_error = (settings["slurm_error"] or "").replace("%x", job_name)
            slurm_settings = {
                "sample_id": settings["sample_id"],
                "stage": stage,
                "slurm_partition": settings["slurm_partition"],
                "slurm_mem": settings["slurm_mem"],
                "slurm_cpus_per_task": settings["slurm_cpus_per_task"],
                "slurm_output": methscan_output,
                "slurm_error": methscan_error,
                "job_name": job_name,
            }
            print(f"[emseq.make_cmd] runner={settings['runner']}")
            print(f"[emseq.make_cmd] stage={stage}")
            print(f"[emseq.make_cmd] sample_id={settings['sample_id']}")
            print(f"[emseq.make_cmd] script={script_path}")
            print(f"[emseq.make_cmd] command={command}")
            if not settings["dry_run"]:
                generate_slurm_script(command, script_path, log_dir, slurm_settings)
                print(f"[emseq.make_cmd] generated={script_path}")
            if settings["submit"] and not settings["dry_run"]:
                submit_script(script_path, settings["runner"])
                print("[emseq.make_cmd] submitted_count=1")
    else:
        raise ValueError(f"unsupported stage for EMSeq entry: {stage}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

