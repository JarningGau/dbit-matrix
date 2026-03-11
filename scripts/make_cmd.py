#!/usr/bin/env python3
"""Generate and optionally submit DBiT workflow commands."""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate executable command scripts for DBiT workflow."
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
        choices=["fastp_split", "demux_extract_bc", "align", "pool", "split", "call"],
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
        help="fastp executable path or command name. Default: fastp.",
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
        help="bwa executable path or command name. Default: bwa.",
    )
    parser.add_argument(
        "--sinto-bin",
        help="sinto executable path or command name. Default: sinto.",
    )
    parser.add_argument(
        "--samtools-bin",
        help="samtools executable path or command name. Default: samtools.",
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
        "--call-caller-script",
        help="Path to methy_caller script for call stage. Default: scripts/methy_caller.py.",
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


def build_demux_local_batch_command(args: argparse.Namespace, sample_work: Path) -> str:
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


def build_call_command(
    args: argparse.Namespace,
    sample_work: Path,
    mode: str,
    spike_name: str | None = None,
) -> str:
    script_path = Path("scripts/call.py")
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
        "--caller-script",
        args.call_caller_script,
    ]
    if args.call_sample_size is not None:
        command.extend(["--sample-size", str(args.call_sample_size)])
    if mode in ("all", "spike"):
        for item in args.spike_in_index:
            command.extend(["--spike-reference", item])
    if spike_name:
        command.extend(["--spike-in-name", spike_name])
    return quoted(command)


def discover_demux_chunks(sample_work: Path) -> list[tuple[str, Path, Path, Path]]:
    chunk_dir = sample_work / "shard_fastq"
    demux_dir = sample_work / "demux"
    chunks: list[tuple[str, Path, Path, Path]] = []
    for r1 in sorted(chunk_dir.glob("*.R1.fq.gz")):
        chunk = r1.name[: -len(".R1.fq.gz")]
        r2 = chunk_dir / f"{chunk}.R2.fq.gz"
        if not r2.exists():
            raise ValueError(f"missing paired R2 for chunk '{chunk}': {r2}")
        out_prefix = demux_dir / chunk
        chunks.append((chunk, r1, r2, out_prefix))
    return chunks


def discover_align_chunks(sample_work: Path) -> list[tuple[str, Path, Path]]:
    demux_dir = sample_work / "demux"
    chunks: list[tuple[str, Path, Path]] = []
    for r1 in sorted(demux_dir.glob("*.R1.demux.fq.gz")):
        chunk = r1.name[: -len(".R1.demux.fq.gz")]
        r2 = demux_dir / f"{chunk}.R2.demux.fq.gz"
        if not r2.exists():
            raise ValueError(f"missing paired R2 for chunk '{chunk}': {r2}")
        chunks.append((chunk, r1, r2))
    return chunks


def discover_call_host_bams(sample_work: Path) -> list[Path]:
    split_dir = sample_work / "split_bams"
    return sorted(split_dir.rglob("*.sorted.bam"))


def discover_call_spike_names(sample_work: Path) -> list[str]:
    pooled_dir = sample_work / "pooled"
    names: list[str] = []
    for bam_path in sorted(pooled_dir.glob("pooled.*.sorted.bam")):
        stem = bam_path.name[: -len(".sorted.bam")]
        if stem == "pooled.byCB":
            continue
        if not stem.startswith("pooled."):
            continue
        spike_name = stem.split(".", 1)[1].strip()
        if spike_name and spike_name not in names:
            names.append(spike_name)
    return names


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

    if any(
        key in slurm_cfg_raw
        for key in ("fastp_split", "demux_extract_bc", "align", "pool", "split", "call")
    ):
        stage_slurm_cfg = slurm_cfg_raw.get(stage, {})
    else:
        # Backward compatibility: legacy flat slurm config.
        stage_slurm_cfg = slurm_cfg_raw
    if stage_slurm_cfg is None:
        stage_slurm_cfg = {}
    if not isinstance(stage_slurm_cfg, dict):
        raise ValueError("selected slurm config must be an object")
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
        "call_caller_script": pick(args.call_caller_script, cfg.get("call_caller_script")),
        "spike_in_index": normalize_spike_in_index(
            pick(args.spike_in_index, cfg.get("spike_in_index"))
        ),
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
        "submit": args.submit,
        "dry_run": args.dry_run,
    }

    required = ["runner", "sample_id"]
    if stage == "fastp_split":
        required.extend(["r1", "r2", "number_of_split_parts"])
    elif stage == "demux_extract_bc":
        required.extend(["barcode1_whitelist", "barcode2_whitelist"])
    elif stage == "align":
        required.extend(["bwa_index"])
    elif stage == "pool":
        pass
    elif stage == "split":
        required.extend(["split_barcodes"])
    elif stage == "call":
        required.extend(["call_reference_file", "call_chromosomes"])
    else:
        raise ValueError(f"unsupported stage: {stage}")

    missing = [key for key in required if settings.get(key) in (None, "")]
    if missing:
        raise ValueError(f"missing required settings: {', '.join(missing)}")

    settings["work_root"] = settings["work_root"] or "work"
    settings["fastp_threads"] = settings["fastp_threads"] or 8
    settings["fastp_bin"] = settings["fastp_bin"] or "fastp"
    settings["linker1"] = settings["linker1"] or "GTGGCCGATGTTTCG"
    settings["linker2"] = (
        settings["linker2"] or "ATCCACGTGCTTGAGAGGCCAGAGCATTCG"
    )
    settings["tn5"] = settings["tn5"] or "CATCGGCGTACGACTAGATGTGTATAAGAGACAG"
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
    settings["bwa_threads"] = (
        int(settings["bwa_threads"]) if settings["bwa_threads"] is not None else 2
    )
    if settings["bwa_threads"] <= 0:
        raise ValueError("bwa_threads must be > 0")
    settings["bwa_bin"] = settings["bwa_bin"] or "bwa"
    settings["sinto_bin"] = settings["sinto_bin"] or "sinto"
    settings["samtools_bin"] = settings["samtools_bin"] or "samtools"
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
    settings["call_caller_script"] = (
        settings["call_caller_script"] or "scripts/methy_caller.py"
    )
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


def main() -> int:
    args = parse_args()
    settings = resolve_settings(args)
    sample_work = Path(settings["work_root"]) / settings["sample_id"]
    command_dir = sample_work / "commands"
    log_dir = sample_work / "logs"
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
                module_line="module load fastp",
            )
            generate_slurm_script(command, script_path, log_dir, slurm_args)
        generated_scripts.append(script_path)
    elif settings["stage"] == "demux_extract_bc":
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
            chunks = discover_demux_chunks(sample_work)
            if not chunks:
                if settings["dry_run"]:
                    print(f"[make_cmd] runner={settings['runner']}")
                    print(f"[make_cmd] stage={settings['stage']}")
                    print(f"[make_cmd] sample_id={settings['sample_id']}")
                    print(
                        f"[make_cmd] no chunks found: {sample_work / 'shard_fastq'}/*.R1.fq.gz"
                    )
                    return 0
                raise ValueError(
                    f"no chunks found under: {sample_work / 'shard_fastq'}/*.R1.fq.gz"
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
            chunks = discover_align_chunks(sample_work)
            if not chunks:
                if settings["dry_run"]:
                    print(f"[make_cmd] runner={settings['runner']}")
                    print(f"[make_cmd] stage={settings['stage']}")
                    print(f"[make_cmd] sample_id={settings['sample_id']}")
                    print(
                        f"[make_cmd] no chunks found: {sample_work / 'demux'}/*.R1.demux.fq.gz"
                    )
                    return 0
                raise ValueError(
                    f"no chunks found under: {sample_work / 'demux'}/*.R1.demux.fq.gz"
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
            split_command = build_split_command(command_args, sample_work)
            sort_command = build_split_sort_command(command_args, sample_work)
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] script={split_script_path}")
            print(f"[make_cmd] command={split_command}")
            print(f"[make_cmd] script={sort_script_path}")
            print(f"[make_cmd] command={sort_command}")
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
            generated_scripts.append(split_script_path)
            generated_scripts.append(sort_script_path)
    else:
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
            call_caller_script=settings["call_caller_script"],
            spike_in_index=settings["spike_in_index"],
        )
        if settings["runner"] == "local":
            script_path = command_dir / "06_call.sh"
            command = build_call_command(command_args, sample_work, "all")
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
            host_bams = discover_call_host_bams(sample_work)
            if not host_bams:
                if settings["dry_run"]:
                    print(f"[make_cmd] runner={settings['runner']}")
                    print(f"[make_cmd] stage={settings['stage']}")
                    print(f"[make_cmd] sample_id={settings['sample_id']}")
                    print(
                        f"[make_cmd] no host spot bams found: "
                        f"{sample_work / 'split_bams'}/**/*.sorted.bam"
                    )
                    return 0
                raise ValueError(
                    f"no host spot bams found under: {sample_work / 'split_bams'}/**/*.sorted.bam"
                )
            host_script_path = command_dir / "06_call_host.sbatch"
            host_command = build_call_command(command_args, sample_work, "host")
            spike_names = discover_call_spike_names(sample_work)
            discovered_from_cfg = parse_spike_names(settings["spike_in_index"])
            if discovered_from_cfg:
                spike_names = [name for name in discovered_from_cfg if name in spike_names]
            print(f"[make_cmd] runner={settings['runner']}")
            print(f"[make_cmd] stage={settings['stage']}")
            print(f"[make_cmd] sample_id={settings['sample_id']}")
            print(f"[make_cmd] host_spot_count={len(host_bams)}")
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
            for spike_name in spike_names:
                spike_script_path = command_dir / f"06_call_spike_{spike_name}.sbatch"
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

    for script_path in generated_scripts:
        print(f"[make_cmd] generated={script_path}")

    if settings["submit"]:
        if settings["stage"] == "split" and settings["runner"] == "slurm":
            split_job_id = submit_slurm_script(generated_scripts[0])
            submit_slurm_script(generated_scripts[1], dependency_job_id=split_job_id)
        else:
            for script_path in generated_scripts:
                submit_script(script_path, settings["runner"])
        print(f"[make_cmd] submitted_count={len(generated_scripts)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
