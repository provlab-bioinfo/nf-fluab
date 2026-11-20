#!/usr/bin/env python3

"""
Summarize SeqKit statistics at the sample level.

Input:
    A tab-delimited SeqKit statistics file containing:
        sample
        platform
        read_type
        file
        num_seqs
        sum_len
        min_len
        avg_len
        max_len

Output:
    A tab-delimited sample-level summary containing:
        sample
        platform
        read_type
        num_reads
        total_bases
        min_len
        avg_len
        max_len

Paired-end samples:
    - Must contain exactly one R1 record and one R2 record.
    - R1 and R2 must have the same number of reads.
    - num_reads is the number of read pairs.
    - total_bases is the sum of R1 and R2 bases.
    - avg_len is calculated across individual reads.

Single-end samples:
    - Must contain exactly one record.
    - num_reads is the number of reads.
    - total_bases is the sum of bases.
    - avg_len is calculated across individual reads.

The script stops with an error if paired-end R1 and R2
read counts are different.
"""

import argparse
import csv
import sys
from collections import defaultdict


def parse_args():
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Summarize SeqKit statistics at the sample level."
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input SeqKit TSV file",
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output sample-level summary TSV file",
    )

    return parser.parse_args()


def format_number(value):
    """Format numeric values for the output."""

    if isinstance(value, float):
        if value.is_integer():
            return str(int(value))

        return f"{value:.1f}"

    return str(value)


def get_read_type_from_filename(filename):
    """
    Determine whether a FASTQ file is R1 or R2.

    Supports both:
        sample_R1_001.fastq.gz
        sample_R2_001.fastq.gz

    and DEACON-style:
        sample.deacon_1.fq.gz
        sample.deacon_2.fq.gz

    Returns:
        "R1", "R2", or None
    """

    if "_R1" in filename:
        return "R1"

    if "_R2" in filename:
        return "R2"

    if ".deacon_1." in filename:
        return "R1"

    if ".deacon_2." in filename:
        return "R2"

    return None


def summarize_paired(sample, rows):
    """
    Summarize a paired-end sample.

    A valid paired-end sample must contain exactly one R1
    and one R2 record.

    R1 and R2 must have identical read counts.

    num_reads:
        Number of read pairs.

    total_bases:
        R1 sum_len + R2 sum_len.

    avg_len:
        Total bases divided by the total number of individual reads.

        total_bases / (num_reads * 2)
    """

    r1_rows = []
    r2_rows = []

    for row in rows:

        # Identify R1/R2 from the filename.
        read = get_read_type_from_filename(row["file"])

        if read == "R1":
            r1_rows.append(row)

        elif read == "R2":
            r2_rows.append(row)

        else:
            raise ValueError(
                f"Sample '{sample}' is marked as paired, but file "
                f"'{row['file']}' cannot be identified as R1 or R2."
            )

    # A paired sample must have exactly one R1.
    if len(r1_rows) != 1:
        raise ValueError(
            f"Sample '{sample}' must have exactly one R1 record, "
            f"but found {len(r1_rows)}."
        )

    # A paired sample must have exactly one R2.
    if len(r2_rows) != 1:
        raise ValueError(
            f"Sample '{sample}' must have exactly one R2 record, "
            f"but found {len(r2_rows)}."
        )

    r1 = r1_rows[0]
    r2 = r2_rows[0]

    r1_reads = int(r1["num_seqs"])
    r2_reads = int(r2["num_seqs"])

    # R1 and R2 must contain exactly the same number of reads.
    # Do NOT use min(R1, R2), because that can hide an invalid pair.
    if r1_reads != r2_reads:
        raise ValueError(
            f"Sample '{sample}' has unequal read counts: "
            f"R1={r1_reads}, R2={r2_reads}."
        )

    # Total bases from both mates.
    total_bases = (
        int(r1["sum_len"])
        + int(r2["sum_len"])
    )

    # Minimum read length across both mates.
    min_len = min(
        int(r1["min_len"]),
        int(r2["min_len"]),
    )

    # Maximum read length across both mates.
    max_len = max(
        int(r1["max_len"]),
        int(r2["max_len"]),
    )

    # There are two individual reads per read pair.
    avg_len = total_bases / (r1_reads * 2)

    return {
        "sample": sample,
        "platform": r1["platform"],
        "read_type": "paired",
        "num_reads": r1_reads,
        "total_bases": total_bases,
        "min_len": min_len,
        "avg_len": avg_len,
        "max_len": max_len,
    }


def summarize_single(sample, rows):
    """
    Summarize a single-end sample.

    A single-end sample must contain exactly one record.
    """

    if len(rows) != 1:
        raise ValueError(
            f"Sample '{sample}' is marked as single-end, "
            f"but found {len(rows)} records."
        )

    row = rows[0]

    num_reads = int(row["num_seqs"])
    total_bases = int(row["sum_len"])

    # For single-end reads, every sequence is one individual read.
    avg_len = total_bases / num_reads

    return {
        "sample": sample,
        "platform": row["platform"],
        "read_type": "single",
        "num_reads": num_reads,
        "total_bases": total_bases,
        "min_len": int(row["min_len"]),
        "avg_len": avg_len,
        "max_len": int(row["max_len"]),
    }


def main():
    """Read SeqKit statistics and create the sample-level summary."""

    args = parse_args()

    # Required columns from the SeqKit statistics file.
    required_columns = {
        "sample",
        "platform",
        "read_type",
        "file",
        "num_seqs",
        "sum_len",
        "min_len",
        "max_len",
    }

    # Group all SeqKit records by sample.
    samples = defaultdict(list)

    with open(args.input, "r", newline="") as infile:

        reader = csv.DictReader(
            infile,
            delimiter="\t",
        )

        if reader.fieldnames is None:
            raise ValueError(
                "Input file is empty or has no header."
            )

        # Check for missing required columns.
        missing = required_columns - set(reader.fieldnames)

        if missing:
            raise ValueError(
                "Input file is missing required columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            samples[row["sample"]].append(row)

    results = []

    # Process each sample independently.
    for sample, rows in samples.items():

        # All records belonging to one sample must have the
        # same read_type.
        read_types = {
            row["read_type"]
            for row in rows
        }

        if len(read_types) != 1:
            raise ValueError(
                f"Sample '{sample}' has inconsistent read_type values: "
                f"{', '.join(sorted(read_types))}."
            )

        read_type = next(iter(read_types))

        # All records belonging to one sample must have the
        # same sequencing platform.
        platforms = {
            row["platform"]
            for row in rows
        }

        if len(platforms) != 1:
            raise ValueError(
                f"Sample '{sample}' has inconsistent platform values: "
                f"{', '.join(sorted(platforms))}."
            )

        # Apply paired-end or single-end logic.
        if read_type == "paired":
            result = summarize_paired(
                sample,
                rows,
            )

        elif read_type == "single":
            result = summarize_single(
                sample,
                rows,
            )

        else:
            raise ValueError(
                f"Sample '{sample}' has unsupported "
                f"read_type '{read_type}'."
            )

        results.append(result)

    # Sort the output by sample name.
    results.sort(
        key=lambda x: x["sample"]
    )

    # Final output columns.
    # The original FASTQ filename is intentionally not included.
    fieldnames = [
        "sample",
        "platform",
        "read_type",
        "num_reads",
        "total_bases",
        "min_len",
        "avg_len",
        "max_len",
    ]

    with open(
        args.output,
        "w",
        newline="",
    ) as outfile:

        writer = csv.DictWriter(
            outfile,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()

        for result in results:

            writer.writerow(
                {
                    "sample": result["sample"],
                    "platform": result["platform"],
                    "read_type": result["read_type"],
                    "num_reads": result["num_reads"],
                    "total_bases": result["total_bases"],
                    "min_len": result["min_len"],
                    "avg_len": format_number(
                        result["avg_len"]
                    ),
                    "max_len": result["max_len"],
                }
            )


if __name__ == "__main__":

    try:
        main()

    except Exception as e:

        # Any validation error causes the script to fail.
        # This prevents an incorrect summary from being generated.
        print(
            f"ERROR: {e}",
            file=sys.stderr,
        )

        sys.exit(1)
