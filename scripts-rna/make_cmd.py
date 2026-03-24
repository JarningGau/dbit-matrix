#!/usr/bin/env python3
"""Generate and optionally submit DBiT-RNA workflow commands."""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
from pathlib import Path

STAGE_SEQUENCE = ["demux_extract_bc", "align"]
STAGE_CHOICES = [*STAGE_SEQUENCE, "all"]
STAGE_REQUIRED_FIELDS = {
    "demux_extract_bc": [
        "r1",
        "r2",
        "barcode1_whitelist",
        "barcode2_whitelist",
        "linker_bc",
        "umi_left",
        "umi_len",
    ],
    "align": ["star_bin", "star_genome_dir", "solo_cb_whitelist"],
}


def quoted(items: list[str]) -> str:
    return " ".join(shlex.quote(i) for i in items)


def load_workflow_config(path: str | None) -> dict:
    if not path:
        return {}
    content = Path(path).read_text(encoding="utf-8")
    data = json.loads(content)
    if not isinstance(data, dict):
        raise ValueError("Workflow config must be a JSON object.")
    return data


def pick(cli_value, cfg_value, default=None):
    if cli_value is not None:
        return cli_value
    if cfg_value is not None:
        return cfg_value
    return default


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate executable command scripts for DBiT-RNA workflow."
    )
    parser.add_argument("--workflow-config", help="JSON config path.")
    parser.add_argument("--runner", choices=["local", "slurm"], help="Runner backend.")
    parser.add_argument("--stage", choices=STAGE_CHOICES, help="Workflow stage.")
    parser.add_argument("--sample-id", help="Sample identifier.")
    parser.add_argument("--r1", help="Input R1 FASTQ(.gz).")
    parser.add_argument("--r2", help="Input R2 FASTQ(.gz).")
    parser.add_argument("--work-root", help="Work root directory.")
    parser.add_argument("--barcode1-whitelist", help="Barcode1 whitelist path.")
    parser.add_argument("--barcode2-whitelist", help="Barcode2 whitelist path.")
    parser.add_argument("--linker-bc", help="BC2-BC1 linker sequence.")
    parser.add_argument("--umi-left", help="Left anchor before UMI.")
    parser.add_argument("--umi-len", type=int, help="UMI length.")
    parser.add_argument(
        "--barcode-hamming-distance",
        type=int,
        help="Max Hamming distance for whitelist fallback.",
    )
    parser.add_argument(
        "--linker-edit-distance",
        type=int,
        help="Max edit distance allowed for linker detection.",
    )
    parser.add_argument("--gzip-level", type=int, help="gzip compress level (0-9).")
    parser.add_argument("--star-bin", help="STAR executable path or command name.")
    parser.add_argument("--star-genome-dir", help="STAR genome index directory.")
    parser.add_argument("--solo-cb-whitelist", help="STARsolo --soloCBwhitelist value.")
    parser.add_argument(
        "--solo-cell-filter",
        help="STARsolo cell filter mode. Default: EmptyDrops_CR.",
    )
    parser.add_argument(
        "--out-tmp-dir",
        help="Optional STAR outTmpDir (recommended on WSL).",
    )
    parser.add_argument("--align-threads", type=int, help="STAR thread count.")
    parser.add_argument(
        "--slurm",
        help="Inline slurm json (optional). Usually from workflow config slurm key.",
    )
    parser.add_argument(
        "--submit", action="store_true", help="Submit generated script immediately."
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="Print intended actions without writing."
    )
    return parser.parse_args()


def resolve_settings(args: argparse.Namespace) -> argparse.Namespace:
    cfg = load_workflow_config(args.workflow_config)
    args.runner = pick(args.runner, cfg.get("runner"))
    args.stage = pick(args.stage, cfg.get("stage"), "demux_extract_bc")
    args.sample_id = pick(args.sample_id, cfg.get("sample_id"))
    args.r1 = pick(args.r1, cfg.get("r1"))
    args.r2 = pick(args.r2, cfg.get("r2"))
    args.work_root = pick(args.work_root, cfg.get("work_root"), "work")
    args.barcode1_whitelist = pick(args.barcode1_whitelist, cfg.get("barcode1_whitelist"))
    args.barcode2_whitelist = pick(args.barcode2_whitelist, cfg.get("barcode2_whitelist"))
    args.linker_bc = pick(args.linker_bc, cfg.get("linker_bc"))
    args.umi_left = pick(args.umi_left, cfg.get("umi_left"))
    args.umi_len = pick(args.umi_len, cfg.get("umi_len"), 10)
    args.barcode_hamming_distance = pick(
        args.barcode_hamming_distance, cfg.get("barcode_hamming_distance"), 1
    )
    args.linker_edit_distance = pick(args.linker_edit_distance, cfg.get("linker_edit_distance"), 1)
    args.gzip_level = pick(args.gzip_level, cfg.get("gzip_level"), 1)
    args.star_bin = pick(args.star_bin, cfg.get("star_bin"), "STAR")
    args.star_genome_dir = pick(args.star_genome_dir, cfg.get("star_genome_dir"))
    args.solo_cb_whitelist = pick(args.solo_cb_whitelist, cfg.get("solo_cb_whitelist"))
    args.solo_cell_filter = pick(
        args.solo_cell_filter, cfg.get("solo_cell_filter"), "EmptyDrops_CR"
    )
    args.out_tmp_dir = pick(args.out_tmp_dir, cfg.get("out_tmp_dir"))
    args.align_threads = pick(args.align_threads, cfg.get("align_threads"), 1)

    slurm_cfg = cfg.get("slurm", {})
    if args.slurm:
        slurm_cfg = json.loads(args.slurm)
    if not isinstance(slurm_cfg, dict):
        raise ValueError("slurm config must be a JSON object.")
    args.slurm = slurm_cfg
    return args


def sample_work_dir(args: argparse.Namespace) -> Path:
    return Path(args.work_root) / args.sample_id


def ensure_required(args: argparse.Namespace, stage: str) -> None:
    missing = []
    for key in ["runner", "sample_id", *STAGE_REQUIRED_FIELDS.get(stage, [])]:
        if getattr(args, key, None) in (None, ""):
            missing.append(key)
    if missing:
        raise ValueError(f"Missing required settings for stage '{stage}': {', '.join(missing)}")


def write_script(path: Path, content: str, dry_run: bool) -> None:
    if dry_run:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    path.chmod(0o755)


def build_demux_command(args: argparse.Namespace, sample_work: Path) -> str:
    return quoted(
        [
            sys.executable,
            "scripts-rna/extract_bc.py",
            "--r1",
            args.r1,
            "--r2",
            args.r2,
            "--output-prefix",
            str(sample_work / "demux" / args.sample_id),
            "--barcode1-whitelist",
            args.barcode1_whitelist,
            "--barcode2-whitelist",
            args.barcode2_whitelist,
            "--linker-bc",
            args.linker_bc,
            "--umi-left",
            args.umi_left,
            "--umi-len",
            str(args.umi_len),
            "--barcode-hamming-distance",
            str(args.barcode_hamming_distance),
            "--linker-edit-distance",
            str(args.linker_edit_distance),
            "--gzip-level",
            str(args.gzip_level),
        ]
    )


def build_align_command(args: argparse.Namespace, sample_work: Path) -> str:
    cmd = [
        sys.executable,
        "scripts-rna/align.py",
        "--work-path",
        str(sample_work),
        "--star-bin",
        args.star_bin,
        "--star-genome-dir",
        args.star_genome_dir,
        "--solo-cb-whitelist",
        args.solo_cb_whitelist,
        "--barcode-whitelist",
        args.barcode1_whitelist,
        "--solo-cell-filter",
        args.solo_cell_filter,
        "--threads",
        str(args.align_threads),
    ]
    if args.out_tmp_dir:
        cmd.extend(["--out-tmp-dir", args.out_tmp_dir])
    return quoted(cmd)

def slurm_header(job_name: str, slurm_cfg: dict) -> str:
    lines = [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name={job_name}",
    ]
    if slurm_cfg.get("partition"):
        lines.append(f"#SBATCH --partition={slurm_cfg['partition']}")
    if slurm_cfg.get("mem"):
        lines.append(f"#SBATCH --mem={slurm_cfg['mem']}")
    if slurm_cfg.get("cpus_per_task"):
        lines.append(f"#SBATCH --cpus-per-task={slurm_cfg['cpus_per_task']}")
    lines.extend(["set -euo pipefail", ""])
    return "\n".join(lines)


def make_local_stage_script(path: Path, command: str, dry_run: bool) -> Path:
    script = "#!/usr/bin/env bash\nset -euo pipefail\n\n" + command + "\n"
    write_script(path, script, dry_run)
    return path


def make_slurm_stage_script(
    path: Path, command: str, job_name: str, slurm_cfg: dict, dry_run: bool
) -> Path:
    script = slurm_header(job_name, slurm_cfg) + command + "\n"
    write_script(path, script, dry_run)
    return path


def generate_single_stage(args: argparse.Namespace, stage: str, dry_run: bool) -> Path:
    ensure_required(args, stage)
    sample_work = sample_work_dir(args)
    commands_dir = sample_work / "commands"
    if stage == "demux_extract_bc":
        command = build_demux_command(args, sample_work)
        if args.runner == "local":
            return make_local_stage_script(commands_dir / "01_demux_extract_bc.sh", command, dry_run)
        return make_slurm_stage_script(
            commands_dir / "01_demux_extract_bc.sbatch",
            command,
            f"{args.sample_id}.rna.demux",
            args.slurm.get("demux_extract_bc", {}),
            dry_run,
        )
    if stage == "align":
        command = build_align_command(args, sample_work)
        if args.runner == "local":
            return make_local_stage_script(commands_dir / "02_align.sh", command, dry_run)
        return make_slurm_stage_script(
            commands_dir / "02_align.sbatch",
            command,
            f"{args.sample_id}.rna.align",
            args.slurm.get("align", {}),
            dry_run,
        )
    raise ValueError(f"Unsupported stage: {stage}")


def generate_all(args: argparse.Namespace, dry_run: bool) -> Path:
    sample_work = sample_work_dir(args)
    commands_dir = sample_work / "commands"
    stage_paths = [generate_single_stage(args, stage, dry_run) for stage in STAGE_SEQUENCE]
    if args.runner == "local":
        body = ["#!/usr/bin/env bash", "set -euo pipefail", ""]
        body.extend(str(p) for p in stage_paths)
        body.append("")
        run_path = commands_dir / "run.sh"
        write_script(run_path, "\n".join(body), dry_run)
        return run_path

    run_lines = slurm_header(f"{args.sample_id}.rna.run", args.slurm.get("align", {})).splitlines()
    run_lines.append(f"jid1=$(sbatch {shlex.quote(str(stage_paths[0]))} | awk '{{print $4}}')")
    run_lines.append(
        f"sbatch --dependency=afterok:$jid1 {shlex.quote(str(stage_paths[1]))}"
    )
    run_lines.append("")
    run_path = commands_dir / "run.sbatch"
    write_script(run_path, "\n".join(run_lines), dry_run)
    return run_path


def submit_script(path: Path) -> None:
    if path.suffix == ".sbatch":
        subprocess.run(["sbatch", str(path)], check=True)
        return
    subprocess.run(["bash", str(path)], check=True)


def main() -> int:
    args = resolve_settings(parse_args())
    if args.stage == "all":
        target = generate_all(args, args.dry_run)
    else:
        target = generate_single_stage(args, args.stage, args.dry_run)
    print(target)
    if args.submit and not args.dry_run:
        submit_script(target)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
