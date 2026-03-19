#!/usr/bin/env python3
"""Generate and optionally submit EMSeq fastp_split commands.

This is an EMSeq-only entrypoint that currently supports a single stage:
`fastp_split`. It intentionally reuses the shared `scripts/fastp_split.py`
implementation while keeping EMSeq orchestration separate from the TAPS
workflow.
"""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
from pathlib import Path


STAGE_SEQUENCE = ["fastp_split"]
STAGE_CHOICES = STAGE_SEQUENCE
STAGE_REQUIRED_FIELDS = {
    "fastp_split": ["r1", "r2", "number_of_split_parts"],
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
        description="Generate executable command scripts for EMSeq fastp_split."
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
        raise ValueError("number_of_split_parts is required for fastp_split")
    settings["number_of_split_parts"] = int(settings["number_of_split_parts"])
    if settings["number_of_split_parts"] <= 0:
        raise ValueError("number_of_split_parts must be > 0")

    settings["fastp_bin"] = normalize_executable_setting(
        settings["fastp_bin"], "fastp"
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
    command: str, output_path: Path, log_dir: Path, settings: dict
) -> None:
    log_dir.mkdir(parents=True, exist_ok=True)
    job_name = f"emseq_fastp_split_{settings['sample_id']}"
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


def main() -> int:
    args = parse_args()
    settings = resolve_settings(args)

    sample_work = Path(settings["work_root"]) / settings["sample_id"]
    command_dir = sample_work / "commands"
    log_dir = sample_work / "logs"

    if settings["stage"] != "fastp_split":
        raise ValueError(f"unsupported stage for EMSeq entry: {settings['stage']}")

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
    print(f"[emseq.make_cmd] stage={settings['stage']}")
    print(f"[emseq.make_cmd] sample_id={settings['sample_id']}")
    print(f"[emseq.make_cmd] script={script_path}")
    print(f"[emseq.make_cmd] command={command}")

    if settings["dry_run"]:
        return 0

    if settings["runner"] == "local":
        generate_local_script(command, script_path)
    else:
        generate_slurm_script(command, script_path, log_dir, settings)

    print(f"[emseq.make_cmd] generated={script_path}")

    if settings["submit"]:
        submit_script(script_path, settings["runner"])
        print("[emseq.make_cmd] submitted_count=1")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

