#!/usr/bin/env python3
"""Methylation caller for DBiT TAPS data."""

from __future__ import annotations

import argparse
from pathlib import Path

import pysam


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Methylation caller for TAPS-seq data")
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
        help="Output path for headerless Bismark-like .cov file.",
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
        help="Only process first N CpG sites per chromosome.",
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
        help="Batch size in bp for grouping CpG positions. Default: 10,000,000.",
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
    cpg_positions: list[tuple[int, int]],
    batch_size: int,
) -> list[list[tuple[int, int]]]:
    if not cpg_positions:
        return []
    batches: list[list[tuple[int, int]]] = []
    current_batch: list[tuple[int, int]] = []
    current_start = cpg_positions[0][0]
    for start_pos, end_pos in cpg_positions:
        if start_pos - current_start <= batch_size:
            current_batch.append((start_pos, end_pos))
        else:
            if current_batch:
                batches.append(current_batch)
            current_batch = [(start_pos, end_pos)]
            current_start = start_pos
    if current_batch:
        batches.append(current_batch)
    return batches


def process_pileup_column(
    pileup_column,
    chromosome: str,
    target_flags: set[int],
) -> dict[str, object]:
    tg_count = 0
    ca_count = 0
    cg_count = 0
    total_coverage = 0

    for pileup_read in pileup_column.pileups:
        if pileup_read.alignment.flag not in target_flags:
            continue
        if pileup_read.is_del or pileup_read.is_refskip:
            continue
        if pileup_read.query_position is None:
            continue
        read_seq = pileup_read.alignment.query_sequence
        base = read_seq[pileup_read.query_position]
        next_pos = pileup_read.query_position + 1
        if next_pos >= len(read_seq):
            continue
        dinucleotide = (base + read_seq[next_pos]).upper()
        if dinucleotide == "TG":
            tg_count += 1
            total_coverage += 1
        elif dinucleotide == "CA":
            ca_count += 1
            total_coverage += 1
        elif dinucleotide == "CG":
            cg_count += 1
            total_coverage += 1

    total_mc = tg_count + ca_count
    total_reads = total_mc + cg_count
    methylation_percent = (total_mc / total_reads) * 100 if total_reads > 0 else 0
    return {
        "chrom": chromosome,
        "pos": pileup_column.pos,
        "TG_counts": tg_count,
        "CA_counts": ca_count,
        "CG_counts": cg_count,
        "methylation_percent": round(methylation_percent, 2),
        "coverage": total_coverage,
    }


def process_cpg_batch(
    inbam,
    chromosome: str,
    cpg_batch: list[tuple[int, int]],
    min_base_quality: int,
    min_mapping_quality: int,
    target_flags: set[int],
    max_depth: int,
) -> list[dict[str, object]]:
    if not cpg_batch:
        return []
    batch_start = min(pos[0] for pos in cpg_batch)
    batch_end = max(pos[1] for pos in cpg_batch)
    target_positions = {pos[0] for pos in cpg_batch}
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
        if pileup_column.pos in target_positions:
            results.append(process_pileup_column(pileup_column, chromosome, target_flags))
    return results


def process_all_cpg_positions(
    inbam,
    chromosome: str,
    cpg_positions: list[tuple[int, int]],
    min_base_quality: int,
    min_mapping_quality: int,
    max_depth: int,
    batch_size: int,
) -> list[dict[str, object]]:
    target_flags = {99, 147, 83, 163}
    results: list[dict[str, object]] = []
    batches = create_position_batches(cpg_positions, batch_size)
    for cpg_batch in batches:
        results.extend(
            process_cpg_batch(
                inbam,
                chromosome,
                cpg_batch,
                min_base_quality,
                min_mapping_quality,
                target_flags,
                max_depth,
            )
        )
    return results


def write_results_to_file(results: list[dict[str, object]], output_file: str) -> None:
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        for result in results:
            chrom = result["chrom"]
            pos = result["pos"]
            methy_percent = result["methylation_percent"]
            mc = result["TG_counts"] + result["CA_counts"]
            c_unmeth = result["CG_counts"]
            handle.write(f"{chrom}\t{pos}\t{pos}\t{methy_percent:.2f}\t{mc}\t{c_unmeth}\n")


def load_reference_and_find_cpg(
    reference_file: str,
    chromosome: str,
    sample_size: int | None,
) -> list[tuple[int, int]]:
    with pysam.FastaFile(reference_file) as fasta:
        sequence = fasta.fetch(chromosome).upper()
    cpg_positions: list[tuple[int, int]] = []
    for start in range(len(sequence) - 1):
        if sequence[start] == "C" and sequence[start + 1] == "G":
            cpg_positions.append((start, start + 1))
    if sample_size and len(cpg_positions) > sample_size:
        cpg_positions = cpg_positions[:sample_size]
    return cpg_positions


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
    verbose: bool,
) -> list[dict[str, object]]:
    chromosomes = parse_chromosome_csv(chromosome)
    validate_chromosomes_in_reference(reference_file, chromosomes)
    inbam = open_bam_file(bam_file)
    try:
        all_results: list[dict[str, object]] = []
        for idx, chrom in enumerate(chromosomes, start=1):
            if verbose:
                print(f"[methy_caller] chromosome={chrom} ({idx}/{len(chromosomes)})")
            cpg_positions = load_reference_and_find_cpg(reference_file, chrom, sample_size)
            chrom_results = process_all_cpg_positions(
                inbam=inbam,
                chromosome=chrom,
                cpg_positions=cpg_positions,
                min_base_quality=min_base_quality,
                min_mapping_quality=min_mapping_quality,
                max_depth=max_depth,
                batch_size=batch_size,
            )
            all_results.extend(chrom_results)
        write_results_to_file(all_results, output_file)
        return all_results
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
    chromosomes = parse_chromosome_csv(args.chromosome)
    validate_chromosomes_in_reference(args.reference_file, chromosomes)

    print(f"[methy_caller] reference={args.reference_file}")
    print(f"[methy_caller] bam={args.bam_file}")
    print(f"[methy_caller] chromosome={args.chromosome}")
    print(f"[methy_caller] output={args.output}")
    if args.sample_size is not None:
        print(f"[methy_caller] sample_size={args.sample_size}")

    if args.dry_run:
        print("[methy_caller] dry_run=1")
        return 0

    results = methylation_caller(
        reference_file=args.reference_file,
        bam_file=args.bam_file,
        output_file=args.output,
        chromosome=args.chromosome,
        min_base_quality=args.min_base_quality,
        min_mapping_quality=args.min_mapping_quality,
        sample_size=args.sample_size,
        max_depth=args.max_depth,
        batch_size=args.batch_size,
        verbose=args.verbose,
    )
    covered_sites = len([result for result in results if int(result["coverage"]) > 0])
    print(f"[methy_caller] done processed_sites={len(results)} covered_sites={covered_sites}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
