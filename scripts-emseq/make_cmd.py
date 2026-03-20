#!/usr/bin/env python3
"""Generate and optionally submit EMSeq workflow commands.

This is an EMSeq-only entrypoint that supports:

- `fastp_split`
- `demux_extract_bc`
- `align`
- `pool`
- `split`

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


STAGE_SEQUENCE = ["fastp_split", "demux_extract_bc", "align", "pool", "split"]
STAGE_CHOICES = STAGE_SEQUENCE
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
    script_path = Path("scripts-emseq/extract_bc.py")
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
            "--linker1",
            args.linker1,
            "--linker2",
            args.linker2,
            "--tn5",
            args.tn5,
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
    py = quoted([sys.executable, "scripts-emseq/extract_bc.py"])
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
        f"--linker1 {shlex.quote(args.linker1)} "
        f"--linker2 {shlex.quote(args.linker2)} "
        f"--tn5 {shlex.quote(args.tn5)} "
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
    if any(key in slurm_cfg_raw for key in STAGE_SEQUENCE):
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

    stage_slurm_cfg = select_stage_slurm_cfg(slurm_cfg_raw, stage)

    split_bams_slurm_cfg = resolve_step_slurm_cfg(
        stage_slurm_cfg, "split_bams", {"split_bams", "sort"}
    )
    split_sort_slurm_cfg = resolve_step_slurm_cfg(
        stage_slurm_cfg, "sort", {"split_bams", "sort"}
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
        "spike_in_index": pick(args.spike_in_index, cfg.get("spike_in_index")),
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
    settings["linker1"] = settings["linker1"]
    settings["linker2"] = settings["linker2"]
    settings["tn5"] = settings["tn5"]
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


def main() -> int:
    args = parse_args()
    settings = resolve_settings(args)
    settings["spike_in_index"] = normalize_spike_in_index(settings["spike_in_index"])

    sample_work = Path(settings["work_root"]) / settings["sample_id"]
    command_dir = sample_work / "commands"
    log_dir = sample_work / "logs"

    stage = settings["stage"]

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
            linker1=settings["linker1"],
            linker2=settings["linker2"],
            tn5=settings["tn5"],
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
    else:
        raise ValueError(f"unsupported stage for EMSeq entry: {stage}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

