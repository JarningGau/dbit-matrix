#!/usr/bin/env python3
"""CH methylation caller for DBiT TAPS data."""

from __future__ import annotations

import argparse
from pathlib import Path

import pysam


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="CH methylation caller for TAPS-seq data")
    parser.add_argument(
        "-f",
        "--reference-file",
        required=True,
        help="Path to reference FASTA file.",
    )
    parser.add_argument(
        "-i",
        "--bam-file",
        required=True,
        help="Path to aligned BAM file.",
    )
    parser.add_argument(
        "-c",
        "--chromosome",
        required=True,
        help=(
            "Chromosome names to analyze. "
            "For multiple chromosomes, use comma-separated values."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output path for headerless CH .cov file with trailing context column.",
    )
    parser.add_argument(
        "--min-base-quality",
        type=int,
        default=30,
        help="Minimum base quality score. Default: 30.",
    )
    parser.add_argument(
        "--min-mapping-quality",
        type=int,
        default=1,
        help="Minimum mapping quality score. Default: 1.",
    )
    parser.add_argument(
        "--sample-size",
        type=int,
        default=None,
        help="Only process first N CH sites per chromosome.",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=250,
        help="Maximum read depth for pileup. Default: 250.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=10_000_000,
        help="Batch size in bp for grouping CH positions. Default: 10,000,000.",
    )
    parser.add_argument(
        "--r1-left-trimming",
        type=int,
        default=0,
        help="R1 left-end trimming in read coordinates. Default: 0.",
    )
    parser.add_argument(
        "--r1-right-trimming",
        type=int,
        default=0,
        help="R1 right-end trimming in read coordinates. Default: 0.",
    )
    parser.add_argument(
        "--r2-left-trimming",
        type=int,
        default=0,
        help="R2 left-end trimming in read coordinates. Default: 0.",
    )
    parser.add_argument(
        "--r2-right-trimming",
        type=int,
        default=0,
        help="R2 right-end trimming in read coordinates. Default: 0.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logs.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print resolved parameters only; do not run calling.",
    )
    return parser.parse_args()


def open_bam_file(bam_file: str):
    return pysam.AlignmentFile(bam_file, "rb")


def create_position_batches(
    ch_positions: list[tuple[int, int, str, str]],
    batch_size: int,
) -> list[list[tuple[int, int, str, str]]]:
    if not ch_positions:
        return []
    batches: list[list[tuple[int, int, str, str]]] = []
    current_batch: list[tuple[int, int, str, str]] = []
    current_start = ch_positions[0][0]
    for start_pos, end_pos, context, strand_mode in ch_positions:
        if start_pos - current_start <= batch_size:
            current_batch.append((start_pos, end_pos, context, strand_mode))
        else:
            if current_batch:
                batches.append(current_batch)
            current_batch = [(start_pos, end_pos, context, strand_mode)]
            current_start = start_pos
    if current_batch:
        batches.append(current_batch)
    return batches


def get_forward_query_index(record: pysam.AlignedSegment, query_pos: int, read_len: int) -> int:
    if record.is_reverse:
        return read_len - query_pos - 1
    return query_pos


def is_trimmed_read_position(
    record: pysam.AlignedSegment,
    query_pos: int,
    read_len: int,
    r1_left_trimming: int,
    r1_right_trimming: int,
    r2_left_trimming: int,
    r2_right_trimming: int,
) -> bool:
    if record.is_read1:
        left_trim = r1_left_trimming
        right_trim = r1_right_trimming
    elif record.is_read2:
        left_trim = r2_left_trimming
        right_trim = r2_right_trimming
    else:
        return False
    if record.is_reverse:
        left_cycle = read_len - query_pos
        right_cycle = query_pos + 1
    else:
        left_cycle = query_pos + 1
        right_cycle = read_len - query_pos
    return left_cycle <= left_trim or right_cycle <= right_trim


def process_pileup_column(
    pileup_column,
    chromosome: str,
    context: str,
    strand: str,
    target_flags: set[int],
    r1_left_trimming: int,
    r1_right_trimming: int,
    r2_left_trimming: int,
    r2_right_trimming: int,
) -> dict[str, object]:
    nth_count = 0
    nch_count = 0
    dan_count = 0
    dgn_count = 0
    total_coverage = 0
    for pileup_read in pileup_column.pileups:
        if pileup_read.alignment.flag not in target_flags:
            continue
        if pileup_read.is_del or pileup_read.is_refskip:
            continue
        if pileup_read.query_position is None:
            continue
        record = pileup_read.alignment
        forward_seq = record.get_forward_sequence()
        if not forward_seq:
            continue
        read_len = len(forward_seq)
        if read_len < 2:
            continue
        query_pos = pileup_read.query_position
        forward_index = get_forward_query_index(record, query_pos, read_len)
        if forward_index < 0 or forward_index >= read_len:
            continue
        if is_trimmed_read_position(
            record=record,
            query_pos=query_pos,
            read_len=read_len,
            r1_left_trimming=r1_left_trimming,
            r1_right_trimming=r1_right_trimming,
            r2_left_trimming=r2_left_trimming,
            r2_right_trimming=r2_right_trimming,
        ):
            continue

        base = forward_seq[forward_index].upper()
        if strand == "+":
            if forward_index + 1 >= read_len:
                continue
            neighbor_query_pos = query_pos - 1 if record.is_reverse else query_pos + 1
            if is_trimmed_read_position(
                record=record,
                query_pos=neighbor_query_pos,
                read_len=read_len,
                r1_left_trimming=r1_left_trimming,
                r1_right_trimming=r1_right_trimming,
                r2_left_trimming=r2_left_trimming,
                r2_right_trimming=r2_right_trimming,
            ):
                continue
            next_base = forward_seq[forward_index + 1].upper()
            if next_base != context[1]:
                continue
            if base == "T":
                nth_count += 1
                total_coverage += 1
            elif base == "C":
                nch_count += 1
                total_coverage += 1
            continue

        if strand == "-":
            if base == "A":
                dan_count += 1
                total_coverage += 1
            elif base == "G":
                dgn_count += 1
                total_coverage += 1

    methylated_count = nth_count + dan_count
    unmethylated_count = nch_count + dgn_count
    total_reads = methylated_count + unmethylated_count
    methylation_percent = (methylated_count / total_reads) * 100 if total_reads > 0 else 0
    return {
        "chrom": chromosome,
        "pos": pileup_column.pos,
        "context": context,
        "strand": strand,
        "methylated_count": methylated_count,
        "unmethylated_count": unmethylated_count,
        "methylation_percent": round(methylation_percent, 2),
        "coverage": total_coverage,
    }


def process_ch_batch(
    inbam,
    chromosome: str,
    ch_batch: list[tuple[int, int, str, str]],
    min_base_quality: int,
    min_mapping_quality: int,
    target_flags: set[int],
    max_depth: int,
    r1_left_trimming: int,
    r1_right_trimming: int,
    r2_left_trimming: int,
    r2_right_trimming: int,
) -> list[dict[str, object]]:
    if not ch_batch:
        return []
    batch_start = min(pos[0] for pos in ch_batch)
    batch_end = max(pos[1] for pos in ch_batch)
    target_positions = {pos[0]: (pos[2], pos[3]) for pos in ch_batch}
    results: list[dict[str, object]] = []
    pileups = inbam.pileup(
        chromosome,
        batch_start,
        batch_end,
        ignore_overlaps=True,
        min_base_quality=min_base_quality,
        stepper="samtools",
        redo_baq=False,
        max_depth=max_depth,
        ignore_orphans=True,
        compute_baq=True,
        min_mapping_quality=min_mapping_quality,
    )
    for pileup_column in pileups:
        target = target_positions.get(pileup_column.pos)
        if target is None:
            continue
        context, strand = target
        results.append(
            process_pileup_column(
                pileup_column,
                chromosome,
                context,
                strand,
                target_flags,
                r1_left_trimming,
                r1_right_trimming,
                r2_left_trimming,
                r2_right_trimming,
            )
        )
    return results


def process_all_ch_positions(
    inbam,
    chromosome: str,
    ch_positions: list[tuple[int, int, str, str]],
    min_base_quality: int,
    min_mapping_quality: int,
    max_depth: int,
    batch_size: int,
    r1_left_trimming: int,
    r1_right_trimming: int,
    r2_left_trimming: int,
    r2_right_trimming: int,
    output_handle,
) -> tuple[int, int]:
    target_flags = {99, 147, 83, 163}
    processed_sites = 0
    covered_sites = 0
    batches = create_position_batches(ch_positions, batch_size)
    for ch_batch in batches:
        batch_results = process_ch_batch(
            inbam,
            chromosome,
            ch_batch,
            min_base_quality,
            min_mapping_quality,
            target_flags,
            max_depth,
            r1_left_trimming,
            r1_right_trimming,
            r2_left_trimming,
            r2_right_trimming,
        )
        processed_sites += len(batch_results)
        for result in batch_results:
            if int(result["coverage"]) <= 0:
                continue
            covered_sites += 1
            output_handle.write(format_result_line(result))
    return processed_sites, covered_sites


def format_result_line(result: dict[str, object]) -> str:
    chrom = result["chrom"]
    pos = result["pos"]
    methy_percent = result["methylation_percent"]
    mc = result["methylated_count"]
    c_unmeth = result["unmethylated_count"]
    context = result["context"]
    strand = result["strand"]
    return (
        f"{chrom}\t{pos}\t{pos}\t{methy_percent:.2f}\t{mc}\t{c_unmeth}\t"
        f"{context}\t{strand}\n"
    )


def load_reference_and_find_ch(
    reference_file: str,
    chromosome: str,
    sample_size: int | None,
) -> list[tuple[int, int, str, str]]:
    with pysam.FastaFile(reference_file) as fasta:
        sequence = fasta.fetch(chromosome).upper()
    reverse_context_map = {
        "T": "CA",
        "G": "CC",
        "A": "CT",
    }
    ch_positions: list[tuple[int, int, str, str]] = []
    for start in range(len(sequence) - 1):
        context = sequence[start : start + 2]
        if context in {"CA", "CC", "CT"}:
            ch_positions.append((start, start + 1, context, "+"))
    for pos in range(1, len(sequence) - 1):
        if sequence[pos] != "G":
            continue
        prev_base = sequence[pos - 1]
        context = reverse_context_map.get(prev_base)
        if context is None:
            continue
        ch_positions.append((pos, pos, context, "-"))
    ch_positions.sort(key=lambda item: item[0])
    if sample_size and len(ch_positions) > sample_size:
        ch_positions = ch_positions[:sample_size]
    return ch_positions


def parse_chromosome_csv(chromosome_csv: str) -> list[str]:
    chromosomes = [chrom.strip() for chrom in chromosome_csv.split(",") if chrom.strip()]
    if not chromosomes:
        raise ValueError("chromosome list is empty")
    return chromosomes


def validate_chromosomes_in_reference(reference_file: str, chromosomes: list[str]) -> None:
    with pysam.FastaFile(reference_file) as fasta:
        available = set(fasta.references)
    missing = [chrom for chrom in chromosomes if chrom not in available]
    if missing:
        missing_str = ",".join(missing)
        raise ValueError(
            f"chromosome(s) not found in reference '{reference_file}': {missing_str}"
        )


def methylation_caller(
    reference_file: str,
    bam_file: str,
    output_file: str,
    chromosome: str,
    min_base_quality: int,
    min_mapping_quality: int,
    sample_size: int | None,
    max_depth: int,
    batch_size: int,
    r1_left_trimming: int,
    r1_right_trimming: int,
    r2_left_trimming: int,
    r2_right_trimming: int,
    verbose: bool,
) -> tuple[int, int]:
    chromosomes = parse_chromosome_csv(chromosome)
    validate_chromosomes_in_reference(reference_file, chromosomes)
    inbam = open_bam_file(bam_file)
    try:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        processed_sites = 0
        covered_sites = 0
        with output_path.open("w", encoding="utf-8") as handle:
            for idx, chrom in enumerate(chromosomes, start=1):
                if verbose:
                    print(f"[methy_caller_CH] chromosome={chrom} ({idx}/{len(chromosomes)})")
                ch_positions = load_reference_and_find_ch(reference_file, chrom, sample_size)
                chrom_processed_sites, chrom_covered_sites = process_all_ch_positions(
                    inbam=inbam,
                    chromosome=chrom,
                    ch_positions=ch_positions,
                    min_base_quality=min_base_quality,
                    min_mapping_quality=min_mapping_quality,
                    max_depth=max_depth,
                    batch_size=batch_size,
                    r1_left_trimming=r1_left_trimming,
                    r1_right_trimming=r1_right_trimming,
                    r2_left_trimming=r2_left_trimming,
                    r2_right_trimming=r2_right_trimming,
                    output_handle=handle,
                )
                processed_sites += chrom_processed_sites
                covered_sites += chrom_covered_sites
        return processed_sites, covered_sites
    finally:
        inbam.close()


def main() -> int:
    args = parse_args()
    if args.min_base_quality < 0:
        raise ValueError("min-base-quality must be >= 0")
    if args.min_mapping_quality < 0:
        raise ValueError("min-mapping-quality must be >= 0")
    if args.max_depth <= 0:
        raise ValueError("max-depth must be > 0")
    if args.batch_size <= 0:
        raise ValueError("batch-size must be > 0")
    if args.sample_size is not None and args.sample_size <= 0:
        raise ValueError("sample-size must be > 0 when provided")
    if args.r1_left_trimming < 0:
        raise ValueError("r1-left-trimming must be >= 0")
    if args.r1_right_trimming < 0:
        raise ValueError("r1-right-trimming must be >= 0")
    if args.r2_left_trimming < 0:
        raise ValueError("r2-left-trimming must be >= 0")
    if args.r2_right_trimming < 0:
        raise ValueError("r2-right-trimming must be >= 0")
    chromosomes = parse_chromosome_csv(args.chromosome)
    validate_chromosomes_in_reference(args.reference_file, chromosomes)

    print(f"[methy_caller_CH] reference={args.reference_file}")
    print(f"[methy_caller_CH] bam={args.bam_file}")
    print(f"[methy_caller_CH] chromosome={args.chromosome}")
    print(f"[methy_caller_CH] output={args.output}")
    if args.sample_size is not None:
        print(f"[methy_caller_CH] sample_size={args.sample_size}")
    print(f"[methy_caller_CH] r1_left_trimming={args.r1_left_trimming}")
    print(f"[methy_caller_CH] r1_right_trimming={args.r1_right_trimming}")
    print(f"[methy_caller_CH] r2_left_trimming={args.r2_left_trimming}")
    print(f"[methy_caller_CH] r2_right_trimming={args.r2_right_trimming}")

    if args.dry_run:
        print("[methy_caller_CH] dry_run=1")
        return 0

    processed_sites, covered_sites = methylation_caller(
        reference_file=args.reference_file,
        bam_file=args.bam_file,
        output_file=args.output,
        chromosome=args.chromosome,
        min_base_quality=args.min_base_quality,
        min_mapping_quality=args.min_mapping_quality,
        sample_size=args.sample_size,
        max_depth=args.max_depth,
        batch_size=args.batch_size,
        r1_left_trimming=args.r1_left_trimming,
        r1_right_trimming=args.r1_right_trimming,
        r2_left_trimming=args.r2_left_trimming,
        r2_right_trimming=args.r2_right_trimming,
        verbose=args.verbose,
    )
    print(
        f"[methy_caller_CH] done processed_sites={processed_sites} "
        f"covered_sites={covered_sites}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
