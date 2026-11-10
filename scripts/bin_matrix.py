#!/usr/bin/env python3
"""Build fixed-bin sparse methylation matrices from host per-spot cov files.

For each requested 5mC context, write a 10x-like directory under
`matrix/<bin_label>/<CONTEXT>/` with:
- features.tsv.gz (genomic bins)
- barcodes.tsv.gz (spots)
- matrix.mtx.gz (multi-value sparse: row col meth cov frac)
"""

from __future__ import annotations

import argparse
import gzip
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Iterator, Optional

VALID_CONTEXTS = ("CG", "CA", "CC", "CT", "CH")
CG_SUFFIX = ".CG.cov"
CH_SUFFIX = ".CH.cov"
CH_CONTEXT_SUFFIXES = {
    "CA": ".CA.cov",
    "CC": ".CC.cov",
    "CT": ".CT.cov",
}
CH_DINUCLEOTIDES = frozenset({"CA", "CC", "CT"})


def format_bin_size_label(bin_size: int) -> str:
    """Compact directory name for a bin width in bp.

    Widths below 1000 stay in bp (`500` -> `500bp`).
    Widths from 1000 bp up to 1 Mb use `k` (`10000` -> `10k`, `100000` -> `100k`).
    Widths of 1 Mb and above use `M` (`1000000` -> `1M`).
    Scales that are not whole units keep a trimmed decimal (`1500` -> `1.5k`).
    """
    if bin_size < 1:
        raise ValueError("bin_size must be >= 1")
    if bin_size < 1000:
        return f"{bin_size}bp"
    if bin_size < 1_000_000:
        return _scaled_size_label(bin_size, 1_000, "k")
    return _scaled_size_label(bin_size, 1_000_000, "M")


def _scaled_size_label(bin_size: int, unit: int, suffix: str) -> str:
    scaled = bin_size / unit
    if scaled.is_integer():
        return f"{int(scaled)}{suffix}"
    text = f"{scaled:.6f}".rstrip("0").rstrip(".")
    return f"{text}{suffix}"


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build fixed genomic-bin sparse matrices (meth/cov/frac) from "
            "coverage/host/**/*.{CG,CH}.cov."
        )
    )
    parser.add_argument(
        "--work-path",
        required=True,
        help="Sample work directory containing coverage/host/.",
    )
    parser.add_argument(
        "--reference",
        required=True,
        help="Host FASTA path; sibling .fai is required.",
    )
    parser.add_argument(
        "--bin-size",
        type=int,
        default=100_000,
        help="Fixed genomic bin width in bp. Default: 100000.",
    )
    parser.add_argument(
        "--contexts",
        default="CG,CA",
        help=(
            "Comma-separated 5mC contexts to build: CG,CA,CC,CT,CH. "
            "Default: CG,CA."
        ),
    )
    parser.add_argument(
        "--chromosomes",
        default=None,
        help=(
            "Optional comma-separated contig allowlist. "
            "Default: all contigs in the FASTA .fai."
        ),
    )
    parser.add_argument(
        "--barcodes",
        default=None,
        help="Optional barcode whitelist file (one barcode per line / first column).",
    )
    parser.add_argument(
        "--out-root",
        default="matrix",
        help="Output root under work-path. Default: matrix.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print resolved paths and counts; do not write outputs.",
    )
    return parser.parse_args(argv)


def parse_contexts(raw: str | Iterable[str]) -> list[str]:
    if isinstance(raw, str):
        parts = [p.strip().upper() for p in raw.split(",") if p.strip()]
    else:
        parts = [str(p).strip().upper() for p in raw if str(p).strip()]
    if not parts:
        raise ValueError("at least one context is required")
    unknown = [c for c in parts if c not in VALID_CONTEXTS]
    if unknown:
        raise ValueError(
            f"unsupported contexts {unknown}; allowed: {', '.join(VALID_CONTEXTS)}"
        )
    # Preserve order, drop duplicates.
    seen: set[str] = set()
    out: list[str] = []
    for ctx in parts:
        if ctx not in seen:
            seen.add(ctx)
            out.append(ctx)
    return out


def parse_chromosome_list(raw: Optional[str]) -> Optional[list[str]]:
    if raw is None or str(raw).strip() == "":
        return None
    return [p.strip() for p in str(raw).split(",") if p.strip()]


def read_fai(fai_path: Path) -> list[tuple[str, int]]:
    if not fai_path.is_file():
        raise FileNotFoundError(f"FASTA index required: {fai_path}")
    contigs: list[tuple[str, int]] = []
    with fai_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                raise ValueError(f"malformed .fai line in {fai_path}: {raw_line!r}")
            contigs.append((fields[0], int(fields[1])))
    if not contigs:
        raise ValueError(f"empty FASTA index: {fai_path}")
    return contigs


def build_bins(
    contigs: list[tuple[str, int]],
    bin_size: int,
    chromosomes: Optional[list[str]],
) -> tuple[list[tuple[str, str, int, int]], dict[str, int]]:
    if bin_size < 1:
        raise ValueError("bin_size must be >= 1")
    allow = set(chromosomes) if chromosomes else None
    bins: list[tuple[str, str, int, int]] = []
    chrom_offsets: dict[str, int] = {}
    for chrom, length in contigs:
        if allow is not None and chrom not in allow:
            continue
        if length < 1:
            continue
        chrom_offsets[chrom] = len(bins)
        start = 0
        while start < length:
            end = min(start + bin_size, length)
            bin_id = f"{chrom}:{start}-{end}"
            bins.append((bin_id, chrom, start, end))
            start = end
    if allow is not None:
        missing = [c for c in chromosomes or [] if c not in chrom_offsets]
        if missing:
            raise ValueError(
                "chromosomes not found in FASTA index: " + ", ".join(missing)
            )
    if not bins:
        raise ValueError("no genomic bins resolved from FASTA index / chromosomes")
    return bins, chrom_offsets


def bin_index_for_site(
    chrom_offsets: dict[str, int],
    chrom: str,
    pos: int,
    bin_size: int,
) -> Optional[int]:
    """Map a 0-based cov site coordinate to a global feature index."""
    if chrom not in chrom_offsets or pos < 0:
        return None
    return chrom_offsets[chrom] + (pos // bin_size)


def spot_id_from_cov_path(cov_path: Path) -> str:
    name = cov_path.name
    for suffix in (CG_SUFFIX, *CH_CONTEXT_SUFFIXES.values(), CH_SUFFIX):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    raise ValueError(
        f"expected *{CG_SUFFIX}, *.CA/CC/CT.cov, or *{CH_SUFFIX}, got: {cov_path}"
    )


def load_barcode_whitelist(path: Path) -> set[str]:
    out: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            out.add(line.split("\t")[0].split(",")[0].strip())
    return out


def adapt_whitelist_to_barcodes(
    whitelist: set[str], barcodes: set[str]
) -> set[str]:
    """Adapt a demux-style short-barcode whitelist to concatenated CB IDs.

    Host cov barcodes are ``<bcB><bcA>``. Workflow ``split_barcodes`` is often
    the shared 8-mer (or other short) demux whitelist. When direct intersection
    is empty and every whitelist entry has the same length ``L`` while observed
    barcodes have length ``2L``, expand to the cartesian product ``bcB+bcA``.
    """
    if not whitelist or not barcodes:
        return whitelist
    direct = barcodes & whitelist
    if direct:
        return whitelist
    wl_lens = {len(x) for x in whitelist if x}
    bc_lens = {len(x) for x in barcodes if x}
    if len(wl_lens) != 1 or len(bc_lens) != 1:
        return whitelist
    wl_len = next(iter(wl_lens))
    bc_len = next(iter(bc_lens))
    if wl_len <= 0 or bc_len != 2 * wl_len:
        return whitelist
    return {left + right for left in whitelist for right in whitelist}


def discover_barcodes(
    host_cov_root: Path,
    contexts: list[str],
    whitelist: Optional[set[str]],
) -> list[str]:
    need_cg = "CG" in contexts
    need_ch = any(c in contexts for c in ("CA", "CC", "CT", "CH"))
    barcodes: set[str] = set()
    if need_cg:
        for path in host_cov_root.rglob(f"*{CG_SUFFIX}"):
            barcodes.add(spot_id_from_cov_path(path))
    if need_ch:
        found_split = False
        for suffix in CH_CONTEXT_SUFFIXES.values():
            for path in host_cov_root.rglob(f"*{suffix}"):
                barcodes.add(spot_id_from_cov_path(path))
                found_split = True
        if not found_split:
            for path in host_cov_root.rglob(f"*{CH_SUFFIX}"):
                barcodes.add(spot_id_from_cov_path(path))
    if whitelist is not None:
        adapted = adapt_whitelist_to_barcodes(whitelist, barcodes)
        barcodes &= adapted
    return sorted(barcodes)


def cov_path_for_barcode(
    host_cov_root: Path, barcode: str, suffix: str
) -> Optional[Path]:
    # Prefer layout coverage/host/<cb[0:3]>/<cb>.*.cov; fall back to rglob.
    preferred = host_cov_root / barcode[:3] / f"{barcode}{suffix}"
    if preferred.is_file():
        return preferred
    matches = list(host_cov_root.rglob(f"{barcode}{suffix}"))
    if not matches:
        return None
    if len(matches) > 1:
        raise ValueError(f"multiple cov files for barcode {barcode}: {matches}")
    return matches[0]


def host_cov_paths_for_context(
    host_cov_root: Path, barcode: str, context: str
) -> list[Path]:
    """Resolve host cov files for a context; CH expands to CA/CC/CT when present."""
    if context == "CG":
        path = cov_path_for_barcode(host_cov_root, barcode, CG_SUFFIX)
        return [path] if path is not None else []
    if context in CH_CONTEXT_SUFFIXES:
        split = cov_path_for_barcode(
            host_cov_root, barcode, CH_CONTEXT_SUFFIXES[context]
        )
        if split is not None:
            return [split]
        legacy = cov_path_for_barcode(host_cov_root, barcode, CH_SUFFIX)
        return [legacy] if legacy is not None else []
    if context == "CH":
        split_paths = [
            cov_path_for_barcode(host_cov_root, barcode, suffix)
            for suffix in CH_CONTEXT_SUFFIXES.values()
        ]
        present = [path for path in split_paths if path is not None]
        if present:
            return present
        legacy = cov_path_for_barcode(host_cov_root, barcode, CH_SUFFIX)
        return [legacy] if legacy is not None else []
    return []


def barcode_has_context_inputs(
    host_cov_root: Path, barcode: str, context: str
) -> bool:
    """True when the barcode has any input cov for this context family."""
    if context == "CG":
        return cov_path_for_barcode(host_cov_root, barcode, CG_SUFFIX) is not None
    if any(
        cov_path_for_barcode(host_cov_root, barcode, suffix) is not None
        for suffix in CH_CONTEXT_SUFFIXES.values()
    ):
        return True
    return cov_path_for_barcode(host_cov_root, barcode, CH_SUFFIX) is not None


def iter_cov_counts(
    cov_path: Path,
    context: str,
) -> Iterator[tuple[str, int, int, int]]:
    """Yield (chrom, pos, meth, unmeth) rows matching context."""
    with cov_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if context == "CG":
                if len(fields) < 6:
                    continue
                chrom = fields[0]
                try:
                    pos = int(fields[1])
                    meth = int(fields[4])
                    unmeth = int(fields[5])
                except ValueError:
                    continue
                yield chrom, pos, meth, unmeth
                continue
            if len(fields) < 7:
                continue
            row_context = fields[6]
            if context == "CH":
                if row_context not in CH_DINUCLEOTIDES:
                    continue
            elif row_context != context:
                continue
            chrom = fields[0]
            try:
                pos = int(fields[1])
                meth = int(fields[4])
                unmeth = int(fields[5])
            except ValueError:
                continue
            yield chrom, pos, meth, unmeth


def accumulate_spot_bins(
    cov_path: Path,
    context: str,
    chrom_offsets: dict[str, int],
    bin_size: int,
    n_bins: int,
) -> dict[int, list[int]]:
    acc: dict[int, list[int]] = defaultdict(lambda: [0, 0])
    for chrom, pos, meth, unmeth in iter_cov_counts(cov_path, context):
        idx = bin_index_for_site(chrom_offsets, chrom, pos, bin_size)
        if idx is None or idx < 0 or idx >= n_bins:
            continue
        pair = acc[idx]
        pair[0] += meth
        pair[1] += unmeth
    return acc


def write_features(path: Path, bins: list[tuple[str, str, int, int]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("bin_id\tchrom\tstart\tend\n")
        for bin_id, chrom, start, end in bins:
            handle.write(f"{bin_id}\t{chrom}\t{start}\t{end}\n")


def write_barcodes(path: Path, barcodes: list[str]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for barcode in barcodes:
            handle.write(f"{barcode}\n")


def write_matrix_mtx_gz(
    path: Path,
    n_features: int,
    n_barcodes: int,
    body_path: Path,
    nnz: int,
) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as out_handle:
        out_handle.write("%%MatrixMarket matrix coordinate real general\n")
        out_handle.write("% columns: row col meth cov frac\n")
        out_handle.write(f"{n_features} {n_barcodes} {nnz}\n")
        with body_path.open("r", encoding="utf-8") as body_handle:
            for line in body_handle:
                out_handle.write(line)


def build_context_matrix(
    *,
    host_cov_root: Path,
    context: str,
    barcodes: list[str],
    bins: list[tuple[str, str, int, int]],
    chrom_offsets: dict[str, int],
    bin_size: int,
    out_dir: Path,
    dry_run: bool,
) -> int:
    missing = [
        bc
        for bc in barcodes
        if not barcode_has_context_inputs(host_cov_root, bc, context)
    ]
    if missing and context == "CG":
        raise FileNotFoundError(
            f"missing {CG_SUFFIX} for barcodes (showing up to 5): {missing[:5]}"
        )
    if missing and context != "CG":
        raise FileNotFoundError(
            f"missing CA/CC/CT or {CH_SUFFIX} for barcodes "
            f"(showing up to 5): {missing[:5]}"
        )

    if dry_run:
        print(
            f"[bin_matrix] dry_run context={context} barcodes={len(barcodes)} "
            f"features={len(bins)} out_dir={out_dir}"
        )
        return 0

    out_dir.mkdir(parents=True, exist_ok=True)
    write_features(out_dir / "features.tsv.gz", bins)
    write_barcodes(out_dir / "barcodes.tsv.gz", barcodes)

    nnz = 0
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", delete=False, suffix=".mtx.body"
    ) as body_tmp:
        body_path = Path(body_tmp.name)
        for col_idx, barcode in enumerate(barcodes, start=1):
            cov_paths = host_cov_paths_for_context(host_cov_root, barcode, context)
            acc: dict[int, list[int]] = defaultdict(lambda: [0, 0])
            for cov_path in cov_paths:
                part = accumulate_spot_bins(
                    cov_path,
                    context,
                    chrom_offsets,
                    bin_size,
                    n_bins=len(bins),
                )
                for row_idx, (meth, unmeth) in part.items():
                    slot = acc[row_idx]
                    slot[0] += meth
                    slot[1] += unmeth
            for row_idx in sorted(acc):
                meth, unmeth = acc[row_idx]
                cov = meth + unmeth
                if cov <= 0:
                    continue
                frac = meth / cov
                # 1-indexed feature row
                body_tmp.write(
                    f"{row_idx + 1} {col_idx} {meth} {cov} {frac:.6g}\n"
                )
                nnz += 1

    try:
        write_matrix_mtx_gz(
            out_dir / "matrix.mtx.gz",
            n_features=len(bins),
            n_barcodes=len(barcodes),
            body_path=body_path,
            nnz=nnz,
        )
    finally:
        body_path.unlink(missing_ok=True)

    print(f"[bin_matrix] wrote context={context} nnz={nnz} out_dir={out_dir}")
    return nnz


def run_bin_matrix(
    *,
    work_path: Path,
    reference: Path,
    bin_size: int,
    contexts: list[str] | str,
    chromosomes: Optional[list[str] | str],
    barcodes: Optional[Path],
    out_root: str,
    dry_run: bool,
) -> int:
    work_path = Path(work_path)
    reference = Path(reference)
    ctx_list = parse_contexts(contexts)
    chrom_list = (
        parse_chromosome_list(chromosomes)
        if isinstance(chromosomes, str) or chromosomes is None
        else list(chromosomes)
    )
    host_cov_root = work_path / "coverage" / "host"
    fai_path = Path(str(reference) + ".fai")
    contigs = read_fai(fai_path)
    bins, chrom_offsets = build_bins(contigs, bin_size, chrom_list)

    whitelist = load_barcode_whitelist(Path(barcodes)) if barcodes else None
    barcode_list = discover_barcodes(host_cov_root, ctx_list, whitelist)
    if not barcode_list:
        raise ValueError(f"no barcodes found under {host_cov_root}")

    out_base = work_path / out_root / format_bin_size_label(bin_size)
    print(f"[bin_matrix] work_path={work_path}")
    print(f"[bin_matrix] reference={reference}")
    print(f"[bin_matrix] fai={fai_path}")
    print(f"[bin_matrix] bin_size={bin_size}")
    print(f"[bin_matrix] contexts={','.join(ctx_list)}")
    print(f"[bin_matrix] n_contigs_used={len(chrom_offsets)}")
    print(f"[bin_matrix] n_features={len(bins)}")
    print(f"[bin_matrix] n_barcodes={len(barcode_list)}")
    print(f"[bin_matrix] out_base={out_base}")
    if dry_run:
        print("[bin_matrix] dry_run=1 skip_write")

    for context in ctx_list:
        build_context_matrix(
            host_cov_root=host_cov_root,
            context=context,
            barcodes=barcode_list,
            bins=bins,
            chrom_offsets=chrom_offsets,
            bin_size=bin_size,
            out_dir=out_base / context,
            dry_run=dry_run,
        )
    print("[bin_matrix] done")
    return 0


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    return run_bin_matrix(
        work_path=Path(args.work_path),
        reference=Path(args.reference),
        bin_size=int(args.bin_size),
        contexts=args.contexts,
        chromosomes=args.chromosomes,
        barcodes=Path(args.barcodes) if args.barcodes else None,
        out_root=args.out_root,
        dry_run=bool(args.dry_run),
    )


if __name__ == "__main__":
    raise SystemExit(main())
