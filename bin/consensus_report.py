#!/usr/bin/env python

import argparse
import csv
import re

"""
#id     length  ATCG    RYSWKMBDHV      N
flua_n_2-seg4|H3N2      1701    1647    0       54
flua_n_2-seg5|H3N2      1497    1497    0       0
flua_n_2-seg7|H3N2      982     982     0       0
flua_n_2-seg8|H3N2      838     838     0       0

id_mapping.tsv:

new_id  original_id
flua_n_2-seg4|H3N2      OR327686
flua_n_2-seg5|H3N2      PV063155
...

Coverage TSV is keyed by the original rname (OR327686, PV063155, ...).
"""


def natural_key(s):
    """Split a string into text/number chunks for human-friendly sorting."""
    return [
        int(part) if part.isdigit() else part
        for part in re.split(r"(\d+)", s)
    ]


def read_stats(stats_file):
    stats = {}

    with open(stats_file, "r") as in_f:
        for line in in_f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) < 5:
                continue

            # "id" column from the consensus stats file
            # (now the renamed ID, e.g. flua_n_2-seg4|H3N2)
            segment = fields[0]

            stats[segment] = {
                "length": fields[1],
                "ATCG": fields[2],
                "RYSWKMBDHV": fields[3],
                "N": fields[4],
            }

    return stats


def read_mapping(mapping_file):
    """Read new_id -> original_id mapping TSV."""

    mapping = {}

    with open(mapping_file, "r", newline="") as in_f:
        reader = csv.DictReader(in_f, delimiter="\t")

        for row in reader:
            mapping[row["new_id"]] = row["original_id"]

    return mapping


def read_coverage(coverage_file):
    # Keyed by rname (== original reference ID from the mapping).
    coverage = {}

    with open(coverage_file, "r") as in_f:
        header = None

        for line in in_f:
            line = line.strip()

            if not line:
                continue

            fields = line.split()

            # Identify the header by the #rname column
            if "#rname" in fields:
                header = [
                    field.lstrip("#")
                    for field in fields
                ]
                continue

            if line.startswith("#"):
                continue

            if header is None:
                continue

            if len(fields) != len(header):
                continue

            row = dict(zip(header, fields))

            rname = row["rname"]

            coverage[rname] = {
                "sample": row["sample"],
                "platform": row["platform"],
                "numreads": row["numreads"],
                "meandepth": row["meandepth"],
            }

    return coverage


def compute_completeness(seq_stats):
    """ATCG / length * 100, formatted to 2 decimals, or 'NA'."""
    try:
        length = int(seq_stats["length"])
        atcg = int(seq_stats["ATCG"])

        if length <= 0:
            return "NA"

        return f"{atcg / length * 100:.2f}"
    except (KeyError, ValueError, TypeError):
        return "NA"


def check_stats_consistency(seq_stats):
    """Return True if ATCG + RYSWKMBDHV + N == length."""
    try:
        return (
            int(seq_stats["ATCG"])
            + int(seq_stats["RYSWKMBDHV"])
            + int(seq_stats["N"])
            == int(seq_stats["length"])
        )
    except (KeyError, ValueError, TypeError):
        return False


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Merge consensus sequence statistics and coverage statistics "
            "via a new_id -> original_id mapping."
        )
    )

    parser.add_argument(
        "-s",
        "--stats",
        required=True,
        help="Consensus stats file (keyed by renamed ID)",
    )

    parser.add_argument(
        "-m",
        "--mapping",
        required=True,
        help="Mapping TSV (new_id -> original_id)",
    )

    parser.add_argument(
        "-c",
        "--coverage",
        required=True,
        help="Coverage TSV file (keyed by original rname)",
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output merged consensus report",
    )

    args = parser.parse_args()

    stats = read_stats(args.stats)
    mapping = read_mapping(args.mapping)
    coverage = read_coverage(args.coverage)

    columns = [
        "sample",
        "platform",
        "segment",
        "length",
        "ATCG",
        "RYSWKMBDHV",
        "N",
        "numreads",
        "completeness",
        "meandepth",
    ]

    rows = []

    for segment in stats:
        seq_stats = stats[segment]

        if not check_stats_consistency(seq_stats):
            raise ValueError(
                f"Inconsistent stats for {segment}: "
                f"ATCG + RYSWKMBDHV + N != length"
            )

        # stats id (renamed) -> original id -> coverage rname
        original_id = mapping.get(segment)

        if original_id is None:
            cov = {}
        else:
            cov = coverage.get(original_id, {})

        values = [
            cov.get("sample", "NA"),
            cov.get("platform", "NA"),
            segment,
            seq_stats["length"],
            seq_stats["ATCG"],
            seq_stats["RYSWKMBDHV"],
            seq_stats["N"],
            cov.get("numreads", "NA"),
            compute_completeness(seq_stats),
            cov.get("meandepth", "NA"),
        ]

        rows.append(values)

    # Sample: alphabetical. Segment: natural (seg2 < seg10, seg4|H1N1 < seg4|H3N2).
    # "NA" samples sort last so rows without coverage don't jump to the top.
    rows.sort(
        key=lambda r: (
            r[0] == "NA",
            r[0],
            natural_key(r[2]),
        )
    )

    with open(args.output, "w") as out_f:
        out_f.write("\t".join(columns) + "\n")

        for values in rows:
            out_f.write("\t".join(values) + "\n")


if __name__ == "__main__":
    main()
