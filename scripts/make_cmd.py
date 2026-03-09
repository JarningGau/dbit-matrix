#!/usr/bin/env python3
"""Generate and optionally submit fastp workflow commands."""

from __future__ import annotations

import argparse
import json
import shlex
import shutil
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
        "--submit",
        action="store_true",
        help="Submit immediately after generating command file.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print command and output path without writing files.",
    )

    # Slurm-only options with pragmatic defaults.
    parser.add_argument("--slurm-partition")
    parser.add_argument("--slurm-time")
    parser.add_argument("--slurm-mem")
    parser.add_argument("--slurm-cpus-per-task", type=int)
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


def resolve_settings(args: argparse.Namespace) -> dict:
    cfg: dict = {}
    if args.workflow_config:
        cfg = load_workflow_config(args.workflow_config)

    slurm_cfg = cfg.get("slurm", {})
    if slurm_cfg is None:
        slurm_cfg = {}
    if not isinstance(slurm_cfg, dict):
        raise ValueError("workflow config key 'slurm' must be an object")

    settings = {
        "runner": pick(args.runner, cfg.get("runner")),
        "sample_id": pick(args.sample_id, cfg.get("sample_id")),
        "r1": pick(args.r1, cfg.get("r1")),
        "r2": pick(args.r2, cfg.get("r2")),
        "work_root": pick(args.work_root, cfg.get("work_root")),
        "fastp_threads": pick(args.fastp_threads, cfg.get("fastp_threads")),
        "number_of_split_parts": pick(
            args.number_of_split_parts, cfg.get("number_of_split_parts")
        ),
        "fastp_bin": pick(args.fastp_bin, cfg.get("fastp_bin")),
        "slurm_partition": pick(args.slurm_partition, slurm_cfg.get("partition")),
        "slurm_time": pick(args.slurm_time, slurm_cfg.get("time")),
        "slurm_mem": pick(args.slurm_mem, slurm_cfg.get("mem")),
        "slurm_cpus_per_task": pick(
            args.slurm_cpus_per_task, slurm_cfg.get("cpus_per_task")
        ),
        "submit": args.submit,
        "dry_run": args.dry_run,
    }

    required = ["runner", "sample_id", "r1", "r2", "number_of_split_parts"]
    missing = [key for key in required if settings.get(key) in (None, "")]
    if missing:
        raise ValueError(f"missing required settings: {', '.join(missing)}")

    settings["work_root"] = settings["work_root"] or "work"
    settings["fastp_threads"] = settings["fastp_threads"] or 8
    settings["fastp_bin"] = settings["fastp_bin"] or "fastp"
    settings["slurm_partition"] = settings["slurm_partition"] or "cpu"
    settings["slurm_time"] = settings["slurm_time"] or "02:00:00"
    settings["slurm_mem"] = settings["slurm_mem"] or "16G"
    settings["slurm_cpus_per_task"] = settings["slurm_cpus_per_task"] or 8

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
    content = (
        "#!/usr/bin/env bash\n"
        f"#SBATCH --job-name=dbit_fastp_{args.sample_id}\n"
        f"#SBATCH --partition={args.slurm_partition}\n"
        f"#SBATCH --time={args.slurm_time}\n"
        f"#SBATCH --cpus-per-task={args.slurm_cpus_per_task}\n"
        f"#SBATCH --mem={args.slurm_mem}\n"
        f"#SBATCH --output={log_dir}/%x_%j.out\n"
        f"#SBATCH --error={log_dir}/%x_%j.err\n\n"
        "set -euo pipefail\n\n"
        "module load fastp\n\n"
        f"{command}\n"
    )
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
    print(f"[make_cmd] sample_id={settings['sample_id']}")
    print(f"[make_cmd] script={script_path}")
    print(f"[make_cmd] command={command}")

    if settings["dry_run"]:
        return 0

    if settings["runner"] == "local":
        generate_local_script(command, script_path)
    else:
        slurm_args = argparse.Namespace(
            sample_id=settings["sample_id"],
            slurm_partition=settings["slurm_partition"],
            slurm_time=settings["slurm_time"],
            slurm_mem=settings["slurm_mem"],
            slurm_cpus_per_task=settings["slurm_cpus_per_task"],
        )
        generate_slurm_script(command, script_path, log_dir, slurm_args)

    print(f"[make_cmd] generated={script_path}")

    if settings["submit"]:
        submit_script(script_path, settings["runner"])
        print("[make_cmd] submitted")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
