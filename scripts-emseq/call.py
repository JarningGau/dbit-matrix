#!/usr/bin/env python3
"""EMSeq call stage: biscuit pileup -> VCF -> TAPS-compatible CG coverage.

Keeps intermediates:
  work/<sample>/pileup/**/*.vcf.gz (+ .tbi)

Produces (TAPS-compatible) outputs:
  work/<sample>/coverage/host/<X_index>/<X_index>_<Y_index>.CG.cov
  work/<sample>/coverage/host_mito.CG.cov
  work/<sample>/coverage/<spike_name>.CG.cov

``host_mito.CG.cov`` prefers ``qc/mbias/host.subsampled.sorted.bam`` (from the
``mbias`` stage) when present; otherwise mitochondrial sites are aggregated from
per-spot host coverage after stripping ``mito_chromosomes`` from those files.
"""

from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from scripts.host_subsample_bam import HOST_SUBSAMPLE_SEED  # type: ignore


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run EMSeq methylation calling (pileup + VCF-to-CG coverage) for "
            "host spots and optional spike-ins."
        )
    )
    parser.add_argument(
        "--work-path",
        required=True,
        help="Sample work directory containing split_bams/, pooled/, and qc/mbias/.",
    )
    parser.add_argument(
        "--mode",
        choices=["all", "host", "spike"],
        default="all",
        help="Run host + spike, or host/spike only. Default: all.",
    )
    parser.add_argument(
        "--reference-file",
        required=True,
        help="Host reference FASTA passed to biscuit pileup/mergecg.",
    )
    parser.add_argument(
        "--spike-reference",
        action="append",
        default=[],
        help="Spike reference in NAME=FASTA format. May be specified multiple times.",
    )
    parser.add_argument(
        "--spike-in-name",
        action="append",
        default=[],
        help=(
            "Spike name to call (may be specified multiple times). "
            "If omitted in spike mode, auto-discover from pooled/."
        ),
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=8,
        help="Maximum per-stage pileup jobs to run concurrently. Default: 8.",
    )
    parser.add_argument(
        "--host-threads",
        type=int,
        default=8,
        help="Threads for biscuit pileup in host mode (bgzip uses 1). Default: 8.",
    )
    parser.add_argument(
        "--spike-threads",
        type=int,
        default=8,
        help="Threads for biscuit pileup in spike mode (bgzip uses 1). Default: 8.",
    )
    parser.add_argument(
        "--biscuit-bin",
        default="biscuit",
        help="biscuit executable path or command name. Default: biscuit.",
    )
    parser.add_argument(
        "--bgzip-bin",
        default="bgzip",
        help="bgzip executable path or command name. Default: bgzip.",
    )
    parser.add_argument(
        "--tabix-bin",
        default="tabix",
        help="tabix executable path or command name. Default: tabix.",
    )
    parser.add_argument(
        "--samtools-bin",
        default="samtools",
        help="samtools executable path or command name. Default: samtools.",
    )
    parser.add_argument(
        "--samtools-threads",
        type=int,
        default=4,
        help="Threads for samtools view/sort/index while preparing host mito BAM. Default: 4.",
    )
    parser.add_argument(
        "--host-subsample-fraction",
        type=float,
        default=0.1,
        help="Host subsampling fraction in (0, 1] for host mito BAM fallback. Default: 0.1.",
    )
    parser.add_argument(
        "--host-subsample-seed",
        type=int,
        default=HOST_SUBSAMPLE_SEED,
        help=f"Host subsampling seed for host mito BAM fallback. Default: {HOST_SUBSAMPLE_SEED}.",
    )
    parser.add_argument(
        "--mito-chromosomes",
        default="chrM",
        help=(
            "Comma-separated host contigs treated as mitochondrial. "
            "They are removed from per-spot coverage/host/**/*.CG.cov; "
            "host_mito.CG.cov is built from qc/mbias/host.subsampled.sorted.bam "
            "when that file exists, else aggregated from stripped mito rows. "
            "Default: chrM."
        ),
    )
    parser.add_argument(
        "--call-left-trimming",
        type=int,
        default=None,
        metavar="N",
        help=(
            "If set, pass -5 N to biscuit pileup (minimum distance to read 5' end). "
            "If omitted, biscuit default applies."
        ),
    )
    parser.add_argument(
        "--call-right-trimming",
        type=int,
        default=None,
        metavar="N",
        help=(
            "If set, pass -3 N to biscuit pileup (minimum distance to read 3' end). "
            "If omitted, biscuit default applies."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands only; do not execute calling.",
    )
    return parser.parse_args()


def parse_spike_references(items: list[str]) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for item in items:
        name, sep, reference_path = item.partition("=")
        if not sep or not name.strip() or not reference_path.strip():
            raise ValueError(f"invalid --spike-reference value (expected NAME=FASTA): {item}")
        parsed[name.strip()] = reference_path.strip()
    return parsed


def parse_chromosome_csv(chromosome_csv: str) -> frozenset[str]:
    chromosomes = frozenset(
        item.strip() for item in chromosome_csv.split(",") if item.strip()
    )
    if not chromosomes:
        raise ValueError("--mito-chromosomes resolved to an empty list")
    return chromosomes


def discover_host_spot_bams(split_dir: Path) -> list[Path]:
    return sorted(split_dir.rglob("*.sorted.bam"))


def discover_spike_bams(pooled_dir: Path) -> dict[str, Path]:
    spikes: dict[str, Path] = {}
    for bam in sorted(pooled_dir.glob("pooled.*.sorted.bam")):
        stem = bam.name[: -len(".sorted.bam")]
        if stem == "pooled.byCB":
            continue
        if not stem.startswith("pooled."):
            continue
        spike_name = stem.split(".", 1)[1]
        if spike_name:
            spikes[spike_name] = bam
    return spikes


def run_command(command: list[str], dry_run: bool) -> None:
    print(f"[emseq.call] command={' '.join(map(str, command))}")
    if dry_run:
        return
    subprocess.run(command, check=True)


def ensure_pileup_vcf(
    *,
    args: argparse.Namespace,
    reference_file: str,
    bam_file: Path,
    out_vcf_gz: Path,
    threads: int,
) -> None:
    out_tbi = Path(str(out_vcf_gz) + ".tbi")
    out_vcf_gz.parent.mkdir(parents=True, exist_ok=True)

    if out_vcf_gz.exists():
        if out_tbi.exists():
            print(f"[emseq.call] skip_existing_pileup={out_vcf_gz}")
            return
        tabix_cmd = [args.tabix_bin, "-f", "-p", "vcf", str(out_vcf_gz)]
        run_command(tabix_cmd, args.dry_run)
        return

    pileup_cmd = [
        args.biscuit_bin,
        "pileup",
        "-m",
        "0",
        "-a",
        "0",
        "-c",
        "-u",
        "-p",
    ]
    if args.call_left_trimming is not None:
        pileup_cmd.extend(["-5", str(args.call_left_trimming)])
    if args.call_right_trimming is not None:
        pileup_cmd.extend(["-3", str(args.call_right_trimming)])
    pileup_cmd.extend(
        [
            "-@",
            str(threads),
            reference_file,
            str(bam_file),
        ]
    )
    bgzip_cmd = [args.bgzip_bin, "-@", "1", "-c"]

    print(f"[emseq.call] pileup_bam={bam_file} -> {out_vcf_gz}")
    if args.dry_run:
        pipe = (
            f"{shlex.join(pileup_cmd)} | {shlex.join(bgzip_cmd)} "
            f"> {shlex.quote(str(out_vcf_gz))}"
        )
        print(f"[emseq.call] command_pipe={pipe}")
        return

    pileup_proc = subprocess.Popen(pileup_cmd, stdout=subprocess.PIPE)
    if pileup_proc.stdout is None:
        raise RuntimeError("failed to capture biscuit pileup stdout")
    try:
        with out_vcf_gz.open("wb") as out_handle:
            bgzip_proc = subprocess.Popen(
                bgzip_cmd,
                stdin=pileup_proc.stdout,
                stdout=out_handle,
            )
            pileup_proc.stdout.close()
            bgzip_rc = bgzip_proc.wait()
            pileup_rc = pileup_proc.wait()
            if bgzip_rc != 0 or pileup_rc != 0:
                raise subprocess.CalledProcessError(
                    bgzip_rc if bgzip_rc != 0 else pileup_rc, bgzip_cmd
                )
    finally:
        try:
            if pileup_proc.poll() is None:
                pileup_proc.kill()
        except Exception:
            pass

    tabix_cmd = [args.tabix_bin, "-f", "-p", "vcf", str(out_vcf_gz)]
    run_command(tabix_cmd, args.dry_run)


def ensure_coverage_from_vcf(
    *,
    args: argparse.Namespace,
    reference_file: str,
    vcf_gz: Path,
    out_cov: Path,
) -> None:
    if out_cov.exists():
        print(f"[emseq.call] skip_existing_coverage={out_cov}")
        return
    out_cov.parent.mkdir(parents=True, exist_ok=True)

    vcf2bed_cmd = [args.biscuit_bin, "vcf2bed", "-t", "cg", str(vcf_gz)]
    mergecg_cmd = [args.biscuit_bin, "mergecg", "-c", reference_file, "-"]

    print(f"[emseq.call] vcf2bed+mergecg {vcf_gz} -> {out_cov}")
    if args.dry_run:
        pipe = (
            f"{shlex.join(vcf2bed_cmd)} | {shlex.join(mergecg_cmd)} | "
            f"awk '{{print $1,$2+1,$3-1,$4,$5,$6}}' > {shlex.quote(str(out_cov))}"
        )
        print(f"[emseq.call] command_pipe={pipe}")
        return

    vcf2bed_proc = subprocess.Popen(
        vcf2bed_cmd, stdout=subprocess.PIPE, text=True, encoding="utf-8"
    )
    if vcf2bed_proc.stdout is None:
        raise RuntimeError("failed to capture vcf2bed stdout")
    mergecg_proc = subprocess.Popen(
        mergecg_cmd,
        stdin=vcf2bed_proc.stdout,
        stdout=subprocess.PIPE,
        text=True,
        encoding="utf-8",
    )
    vcf2bed_proc.stdout.close()
    assert mergecg_proc.stdout is not None

    with out_cov.open("w", encoding="utf-8") as out_handle:
        for line in mergecg_proc.stdout:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                continue
            try:
                start = int(fields[1]) + 1
                end = int(fields[2]) - 1
            except ValueError:
                continue
            out_handle.write(
                "\t".join(
                    [fields[0], str(start), str(end), fields[3], fields[4], fields[5]]
                )
                + "\n"
            )

    mergecg_rc = mergecg_proc.wait()
    vcf2bed_rc = vcf2bed_proc.wait()
    if mergecg_rc != 0 or vcf2bed_rc != 0:
        raise subprocess.CalledProcessError(
            mergecg_rc if mergecg_rc != 0 else vcf2bed_rc, mergecg_cmd
        )


def format_methylation_percent(methylated: int, unmethylated: int) -> str:
    total = methylated + unmethylated
    if total <= 0:
        return "0"
    return str(int(round((methylated / total) * 100)))


def strip_mito_from_host_spot_covs_and_collect(
    cov_root: Path,
    mito_chromosomes: frozenset[str],
    *,
    dry_run: bool,
) -> dict[tuple[str, int, int], list[int]]:
    """Remove mitochondrial rows from per-spot host ``*.CG.cov`` files and aggregate them."""
    host_covs = sorted(cov_root.rglob("*.CG.cov"))
    if not host_covs:
        if dry_run:
            print(
                f"[emseq.call] dry_run_skip strip_mito_from_host_spot_covs (no *.CG.cov under "
                f"{cov_root})"
            )
            return {}
        raise ValueError(f"no host coverage files found under: {cov_root}/**/*.CG.cov")

    print(
        f"[emseq.call] strip_mito_from_host_spot_covs count={len(host_covs)} "
        f"mito={sorted(mito_chromosomes)}"
    )
    if dry_run:
        print(
            f"[emseq.call] dry_run_skip in-place strip of mito contigs from per-spot cov files"
        )
        return {}

    mito_rows: dict[tuple[str, int, int], list[int]] = {}
    for cov_path in host_covs:
        kept_lines: list[str] = []
        with cov_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6:
                    continue
                chrom = fields[0]
                if chrom in mito_chromosomes:
                    try:
                        start = int(fields[1])
                        end = int(fields[2])
                        methylated = int(fields[4])
                        unmethylated = int(fields[5])
                    except ValueError:
                        continue
                    key = (chrom, start, end)
                    bucket = mito_rows.setdefault(key, [0, 0])
                    bucket[0] += methylated
                    bucket[1] += unmethylated
                    continue
                kept_lines.append(line)

        with cov_path.open("w", encoding="utf-8") as handle:
            handle.writelines(kept_lines)

    return mito_rows


def write_host_mito_aggregated(
    mito_rows: dict[tuple[str, int, int], list[int]],
    mito_out_cov: Path,
) -> None:
    mito_out_cov.parent.mkdir(parents=True, exist_ok=True)
    with mito_out_cov.open("w", encoding="utf-8") as handle:
        for chrom, start, end in sorted(mito_rows):
            methylated, unmethylated = mito_rows[(chrom, start, end)]
            pct = format_methylation_percent(methylated, unmethylated)
            handle.write(
                "\t".join(
                    [
                        chrom,
                        str(start),
                        str(end),
                        pct,
                        str(methylated),
                        str(unmethylated),
                    ]
                )
                + "\n"
            )


def filter_cov_file_to_mito_chromosomes(
    src_cov: Path,
    dst_cov: Path,
    mito_chromosomes: frozenset[str],
    *,
    dry_run: bool,
) -> None:
    """Keep only rows whose chrom is in ``mito_chromosomes`` (6+ column .CG.cov)."""
    if dry_run:
        print(
            f"[emseq.call] dry_run_skip filter_cov_file_to_mito_chromosomes "
            f"{src_cov} -> {dst_cov}"
        )
        return
    if not src_cov.exists():
        raise ValueError(f"missing intermediate coverage for mito filter: {src_cov}")
    dst_cov.parent.mkdir(parents=True, exist_ok=True)
    with src_cov.open("r", encoding="utf-8") as src_f, dst_cov.open(
        "w", encoding="utf-8"
    ) as dst_f:
        for line in src_f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            if fields[0] in mito_chromosomes:
                dst_f.write(line)


def run_host_mito_from_mbias_bam(
    args: argparse.Namespace,
    work_path: Path,
    mito_chromosomes: frozenset[str],
    mito_out_cov: Path,
) -> None:
    """Pileup ``qc/mbias/host.subsampled.sorted.bam`` and write aggregated mito coverage."""
    mbias_bam = work_path / "qc" / "mbias" / "host.subsampled.sorted.bam"
    pileup_root = work_path / "pileup"
    tmp_cov = work_path / "coverage" / "_host_mito_from_mbias_pooled.CG.cov"
    out_vcf_gz = pileup_root / "host_mito_from_mbias_pooled.vcf.gz"

    print(f"[emseq.call] host_mito_source={mbias_bam}")
    ensure_pileup_vcf(
        args=args,
        reference_file=args.reference_file,
        bam_file=mbias_bam,
        out_vcf_gz=out_vcf_gz,
        threads=args.host_threads,
    )
    ensure_coverage_from_vcf(
        args=args,
        reference_file=args.reference_file,
        vcf_gz=out_vcf_gz,
        out_cov=tmp_cov,
    )
    filter_cov_file_to_mito_chromosomes(
        tmp_cov,
        mito_out_cov,
        mito_chromosomes,
        dry_run=args.dry_run,
    )
    if not args.dry_run and tmp_cov.exists():
        tmp_cov.unlink()


def run_host_spots(args: argparse.Namespace, work_path: Path) -> Path:
    split_dir = work_path / "split_bams"
    host_bams = discover_host_spot_bams(split_dir)
    if not host_bams:
        raise ValueError(f"no host spot BAMs found under: {split_dir}/**/*.sorted.bam")

    print(f"[emseq.call] host_spot_count={len(host_bams)}")
    pileup_root = work_path / "pileup"
    cov_root = work_path / "coverage" / "host"

    def handle_one(bam_file: Path) -> None:
        relative = bam_file.relative_to(split_dir)
        bam_base = bam_file.name[: -len(".sorted.bam")]
        out_vcf_gz = pileup_root / relative.parent / f"{bam_base}.vcf.gz"
        out_cov = cov_root / relative.parent / f"{bam_base}.CG.cov"

        ensure_pileup_vcf(
            args=args,
            reference_file=args.reference_file,
            bam_file=bam_file,
            out_vcf_gz=out_vcf_gz,
            threads=args.host_threads,
        )
        ensure_coverage_from_vcf(
            args=args,
            reference_file=args.reference_file,
            vcf_gz=out_vcf_gz,
            out_cov=out_cov,
        )

    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {executor.submit(handle_one, bam): bam for bam in host_bams}
        for future in as_completed(futures):
            bam = futures[future]
            future.result()
            print(f"[emseq.call] host_bam_done={bam}")
    return cov_root


def run_host_mito(
    args: argparse.Namespace,
    work_path: Path,
    host_cov_root: Path,
    mito_chromosomes: frozenset[str],
) -> None:
    mito_out_cov = work_path / "coverage" / "host_mito.CG.cov"
    mbias_bam = work_path / "qc" / "mbias" / "host.subsampled.sorted.bam"
    host_covs = sorted(host_cov_root.rglob("*.CG.cov"))

    if host_covs:
        mito_rows = strip_mito_from_host_spot_covs_and_collect(
            host_cov_root,
            mito_chromosomes,
            dry_run=args.dry_run,
        )
    else:
        mito_rows = {}
        if args.dry_run:
            print(
                f"[emseq.call] dry_run: no per-spot *.CG.cov under {host_cov_root} yet"
            )

    if mbias_bam.exists():
        run_host_mito_from_mbias_bam(
            args, work_path, mito_chromosomes, mito_out_cov
        )
        print(f"[emseq.call] output_host_mito={mito_out_cov}")
        return

    if args.dry_run:
        print(
            f"[emseq.call] dry_run: would aggregate host_mito from per-spot mito rows "
            f"-> {mito_out_cov} (no {mbias_bam})"
        )
        return

    if not mito_rows and not host_covs:
        raise ValueError(
            f"cannot build host_mito.CG.cov: no per-spot coverage under {host_cov_root} "
            f"and missing mbias host BAM {mbias_bam}"
        )

    write_host_mito_aggregated(mito_rows, mito_out_cov)
    print(f"[emseq.call] output_host_mito={mito_out_cov}")


def run_spikes(args: argparse.Namespace, work_path: Path) -> None:
    pooled_dir = work_path / "pooled"
    discovered_spikes = discover_spike_bams(pooled_dir)
    if not discovered_spikes:
        raise ValueError(
            f"no spike-in BAMs found under: {pooled_dir}/pooled.*.sorted.bam"
        )

    spike_reference_map = parse_spike_references(args.spike_reference)
    requested_names = [n.strip() for n in args.spike_in_name if n.strip()]
    spike_names = requested_names if requested_names else sorted(discovered_spikes.keys())

    print(f"[emseq.call] spike_in_count={len(spike_names)}")
    pileup_root = work_path / "pileup"
    cov_root = work_path / "coverage"

    def handle_one(spike_name: str) -> None:
        spike_bam = discovered_spikes.get(spike_name)
        if spike_bam is None:
            raise ValueError(f"requested spike-in BAM missing for '{spike_name}'")
        spike_reference = spike_reference_map.get(spike_name)
        if spike_reference is None:
            raise ValueError(
                f"missing spike reference for '{spike_name}'. Provide --spike-reference NAME=FASTA."
            )

        bam_base = spike_bam.name[: -len(".sorted.bam")]
        out_vcf_gz = pileup_root / f"{bam_base}.vcf.gz"
        out_cov = cov_root / f"{spike_name}.CG.cov"

        ensure_pileup_vcf(
            args=args,
            reference_file=spike_reference,
            bam_file=spike_bam,
            out_vcf_gz=out_vcf_gz,
            threads=args.spike_threads,
        )
        ensure_coverage_from_vcf(
            args=args,
            reference_file=spike_reference,
            vcf_gz=out_vcf_gz,
            out_cov=out_cov,
        )

    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {executor.submit(handle_one, s): s for s in spike_names}
        for future in as_completed(futures):
            s = futures[future]
            future.result()
            print(f"[emseq.call] spike_done={s}")


def main() -> int:
    args = parse_args()
    if args.jobs <= 0:
        raise ValueError("--jobs must be > 0")
    if args.host_threads <= 0 or args.spike_threads <= 0:
        raise ValueError("--host-threads/--spike-threads must be > 0")
    if args.samtools_threads <= 0:
        raise ValueError("--samtools-threads must be > 0")
    if args.host_subsample_seed < 0:
        raise ValueError("--host-subsample-seed must be >= 0")
    if args.host_subsample_fraction <= 0 or args.host_subsample_fraction > 1:
        raise ValueError("--host-subsample-fraction must be in (0, 1]")
    if args.call_left_trimming is not None and args.call_left_trimming < 0:
        raise ValueError("--call-left-trimming must be >= 0")
    if args.call_right_trimming is not None and args.call_right_trimming < 0:
        raise ValueError("--call-right-trimming must be >= 0")

    work_path = Path(args.work_path)
    mito_chromosomes = parse_chromosome_csv(args.mito_chromosomes)
    print(f"[emseq.call] mode={args.mode}")
    print(f"[emseq.call] work_path={work_path}")
    print(f"[emseq.call] mito_chromosomes={sorted(mito_chromosomes)}")

    spike_reference_items = args.spike_reference or []
    effective_mode = args.mode
    if effective_mode == "all" and not spike_reference_items:
        print("[emseq.call] spike_reference empty: run host only")
        effective_mode = "host"

    if effective_mode in ("all", "host"):
        host_cov_root = run_host_spots(args, work_path)
        run_host_mito(args, work_path, host_cov_root, mito_chromosomes)

    if effective_mode in ("all", "spike"):
        if not spike_reference_items:
            raise ValueError("spike mode requires --spike-reference NAME=FASTA")
        run_spikes(args, work_path)

    print("[emseq.call] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
