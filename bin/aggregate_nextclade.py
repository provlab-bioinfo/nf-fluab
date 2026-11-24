#!/usr/bin/env python3

"""
aggregate_nextclade.py

Merge Nextclade TSV results for influenza HA (segment 4) and
NA (segment 6) into two long-format tables plus one wide
per-(sample, platform, subtype) Nextclade summary table.

Input:
    Directory containing Nextclade TSV output files:
        *.nextclade.tsv

    Only TSV is supported. CSV and JSON are not read.

Only Nextclade datasets for:
    Segment 4 -> HA
    Segment 6 -> NA

are included.

Outputs:
    <prefix>-nextclade_ha.long.tsv
        Segment 4 rows, union header of HA files.
        Includes sample and platform columns.

    <prefix>-nextclade_na.long.tsv
        Segment 6 rows, union header of NA files.
        Includes sample and platform columns.

    <prefix>-nextclade_ha_na.short.tsv
        One row per (sample, platform, subtype), wide HA/NA summary.

Long-format output (per gene):
    - One row per surviving input row.
    - Header is source_file, dataset, sample, platform,
      followed by the union of all columns across that gene's
      input files.
    - sample and platform are parsed from seqName.
    - A file that lacks a given column contributes "NA" for
      that column in its rows.
    - Multiple hits per sample are preserved.
    - No deduplication or collapsing across files.

Summary output:
    One row per (sample, platform, subtype).

    Expected seqName format:
        <sample>-seg<4|6>|<platform>|<subtype>

    Examples:
        flua_i_1-seg4|illumina|H3N2
        flua_i_1-seg6|illumina|H3N2
        flua_i_1-seg4|nanopore|H1N1

    Each summary row carries the HA and NA hits matching
    that sample, platform, and subtype.

    Co-infections:
        A sample evaluated against more than one subtype
        produces one row per subtype, separately for each
        platform.

    Missing values:
        Missing gene hits or empty input fields are written
        as "NA".

    Collapsing:
        Multiple Nextclade results for the same seqName are
        collapsed when they have the same classification.
        The first-seen classification wins.

Filtering:
    - Rows with a non-empty "errors" field are excluded.
    - Rows with qc.overallStatus == "bad" are excluded.
    - No other QC filtering.
"""


import argparse
import csv
import glob
import os
import re


GENE_NAMES = {
    "pb2": "PB2",
    "pb1": "PB1",
    "pa": "PA",
    "ha": "HA",
    "np": "NP",
    "na": "NA",
    "mp": "MP",
    "ns": "NS",
}


LONG_FILENAME_TEMPLATE = {
    "HA": "{prefix}-nextclade_ha.long.tsv",
    "NA": "{prefix}-nextclade_na.long.tsv",
}


SUMMARY_FILENAME_TEMPLATE = (
    "{prefix}-nextclade_ha_na.short.tsv"
)


# Explicit columns added to both long-format outputs.
# sample and platform are extracted from seqName.
""" PREFIX_COLUMNS = [
    "source_file",
    "dataset",
    "sample",
    "platform",
] """
PREFIX_COLUMNS = [
    "sample",
    "platform",
    "source_file",
    "dataset",
]

CLASSIFICATION_COLUMNS = [
    "clade",
    "proposedSubclade",
    "subclade",
]


QC_COLUMN = "qc.overallStatus"
BAD_QC_VALUE = "bad"


# Expected seqName:
#   <sample>-seg<segment>|<platform>|<subtype>
#
# Example:
#   sample8-seg4|nanopore|H1N1
#
# The optional trailing text is ignored.
SEQNAME_RE = re.compile(
    r"^(?P<sample>.+?)-seg(?P<seg>[0-9]+)\|"
    r"(?P<platform>[^|\s]+)\|"
    r"(?P<subtype>[^|\s]+)"
)


SUMMARY_GENE_COLUMNS = {
    "HA": {
        "clade": "HA_clade",
        "proposedSubclade": "HA_proposedSubclade",
        "subclade": "HA_subclade",
        QC_COLUMN: "HA_nextclade_qc.overallStatus",
    },
    "NA": {
        "clade": "NA_clade",
        "proposedSubclade": "NA_proposedSubclade",
        "subclade": "NA_subclade",
        QC_COLUMN: "NA_nextclade_qc.overallStatus",
    },
}


SUMMARY_COLUMNS = [
    "sample",
    "platform",
    "ref_subtype",
    "seg4_id",
    "HA_clade",
    "HA_proposedSubclade",
    "HA_subclade",
    "HA_nextclade_qc.overallStatus",
    "seg6_id",
    "NA_clade",
    "NA_proposedSubclade",
    "NA_subclade",
    "NA_nextclade_qc.overallStatus",
]


def parse_seqname(seqname):
    """
    Extract sample, segment, platform, and subtype from seqName.

    Returns:
        dict with sample, seg, platform, subtype
        or None if seqName does not match the expected format.
    """

    match = SEQNAME_RE.match(seqname or "")

    if not match:
        return None

    return match.groupdict()


def is_bad_qc(value):
    """Return True if QC status is bad, ignoring case and whitespace."""

    return (value or "").strip().lower() == BAD_QC_VALUE


def get_dataset_name(filename):
    """
    Extract the Nextclade dataset name from the output filename.

    Example:
        all_consensus.nextclade.nextstrain---flu---h3n2---ha---CY163680.nextclade.tsv

    Returns:
        nextstrain/flu/h3n2/ha/CY163680
    """

    basename = os.path.basename(filename)

    dataset = basename.removeprefix("all_consensus.nextclade.")
    dataset = dataset.removesuffix(".nextclade.tsv")

    dataset = dataset.replace("---", "/")

    return dataset


def get_gene_from_dataset(dataset):
    """Determine influenza gene from the Nextclade dataset name."""

    parts = dataset.lower().split("/")

    for part in parts:
        if part in GENE_NAMES:
            return GENE_NAMES[part]

    return None


def read_headers(filename):
    """Return fieldnames from a Nextclade TSV file."""

    try:
        with open(filename, "r", newline="") as in_fh:
            reader = csv.reader(in_fh, delimiter="\t")
            header = next(reader, None)

    except OSError as exc:
        print(
            f"WARNING: Cannot read {filename}: {exc}",
            flush=True,
        )
        return []

    return header or []


def collect_union_columns(files):
    """
    Return the union of data columns across files,
    preserving first-seen order.

    Columns that collide with PREFIX_COLUMNS are excluded.
    """

    prefix_set = set(PREFIX_COLUMNS)
    columns = []
    seen = set()

    for filename, _gene in files:
        for field in read_headers(filename):
            if not field:
                continue

            if field in prefix_set:
                continue

            if field in seen:
                continue

            seen.add(field)
            columns.append(field)

    return columns


def classification_key(row):
    """Return the classification tuple used for summary collapsing."""

    return tuple(
        (row.get(column) or "").strip()
        for column in CLASSIFICATION_COLUMNS
    )


def aggregate_gene(files, union_columns):
    """
    Build one output row per surviving input row for a single gene.

    sample and platform are extracted from seqName.

    No collapsing happens here. Multiple hits are preserved.
    """

    rows = []

    total_rows = 0
    kept_rows = 0
    skipped_errors = 0
    skipped_bad_qc = 0
    skipped_unparseable = 0

    for filename, _gene in files:
        dataset = get_dataset_name(filename)

        print(f"Reading: {filename}", flush=True)

        with open(filename, "r", newline="") as in_fh:
            reader = csv.DictReader(in_fh, delimiter="\t")

            if not reader.fieldnames:
                print(
                    f"WARNING: No header found: {filename}",
                    flush=True,
                )
                continue

            if "seqName" not in reader.fieldnames:
                print(
                    f"WARNING: No 'seqName' column found: "
                    f"{filename}; skipping file",
                    flush=True,
                )
                continue

            if "errors" not in reader.fieldnames:
                print(
                    f"WARNING: No 'errors' column found: "
                    f"{filename}; error filtering disabled",
                    flush=True,
                )

            if QC_COLUMN not in reader.fieldnames:
                print(
                    f"WARNING: No '{QC_COLUMN}' column found: "
                    f"{filename}; bad-QC filtering disabled",
                    flush=True,
                )

            for row in reader:
                total_rows += 1

                # 1. Exclude rows with Nextclade errors.
                errors = (row.get("errors") or "").strip()

                if errors:
                    skipped_errors += 1
                    continue

                # 2. Exclude rows with bad overall QC.
                if is_bad_qc(row.get(QC_COLUMN)):
                    skipped_bad_qc += 1
                    continue

                # 3. Extract sample and platform from seqName.
                seqname = row.get("seqName", "")
                parsed = parse_seqname(seqname)

                if parsed is None:
                    skipped_unparseable += 1
                    print(
                        f"WARNING: Cannot parse seqName: "
                        f"{seqname!r} in {filename}",
                        flush=True,
                    )
                    continue

                out = {
                    "source_file": os.path.basename(filename),
                    "dataset": dataset,
                    "sample": parsed["sample"],
                    "platform": parsed["platform"],
                }

                for column in union_columns:
                    value = (row.get(column) or "").strip()
                    out[column] = value if value else "NA"

                rows.append(out)
                kept_rows += 1

    print(
        f"  rows read={total_rows}, kept={kept_rows}, "
        f"with errors={skipped_errors}, "
        f"bad QC={skipped_bad_qc}, "
        f"unparseable seqName={skipped_unparseable}",
        flush=True,
    )

    return rows


def write_gene(rows, union_columns, output_file):
    """
    Write one gene's long-format table.

    Column order:
        source_file, dataset, sample, platform, union columns
    """

    output_dir = os.path.dirname(output_file)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    fieldnames = PREFIX_COLUMNS + union_columns

    rows_sorted = sorted(
        rows,
        key=lambda row: (
            row["sample"],
            row["platform"],
            row["source_file"],
            row.get("seqName", ""),
        ),
    )

    with open(output_file, "w", newline="") as out_fh:
        writer = csv.DictWriter(
            out_fh,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(rows_sorted)

    print(
        f"Wrote {len(rows_sorted)} rows, {len(fieldnames)} columns "
        f"-> {output_file}",
        flush=True,
    )


def write_summary(ha_rows, na_rows, output_file):
    """
    Build the wide summary.

    One row per (sample, platform, subtype).

    HA and NA hits are matched using all three fields:
        sample, platform, subtype

    Multiple results for the same seqName and classification
    are collapsed. The first-seen classification wins.
    """

    # Per-gene collapse:
    # (seqName, classification) -> first row
    collapsed = {
        "HA": {},
        "NA": {},
    }

    counts = {
        "HA": 0,
        "NA": 0,
    }

    for gene, rows in (("HA", ha_rows), ("NA", na_rows)):
        for row in rows:
            seqname = row.get("seqName", "")

            key = (
                seqname,
                classification_key(row),
            )

            if key in collapsed[gene]:
                counts[gene] += 1
                continue

            collapsed[gene][key] = row

        if counts[gene]:
            print(
                f"  summary {gene}: collapsed "
                f"{counts[gene]} duplicate "
                f"(seqName, classification) rows from "
                f"{len(rows)} long-format rows",
                flush=True,
            )

    # Collect hits by (sample, platform, subtype).
    by_key = {}
    order = []

    for gene in ("HA", "NA"):
        expected_seg = "4" if gene == "HA" else "6"

        for row in collapsed[gene].values():
            seqname = row.get("seqName", "")
            parsed = parse_seqname(seqname)

            if parsed is None:
                print(
                    f"WARNING: Cannot parse seqName for summary: "
                    f"{seqname!r}",
                    flush=True,
                )
                continue

            if parsed["seg"] != expected_seg:
                print(
                    f"WARNING: Segment mismatch for "
                    f"{parsed['sample']}: got seg{parsed['seg']} "
                    f"in {gene} file; skipping {seqname!r}",
                    flush=True,
                )
                continue

            sample = parsed["sample"]
            platform = parsed["platform"]
            subtype = parsed["subtype"]

            key = (sample, platform, subtype)

            entry = by_key.get(key)

            if entry is None:
                entry = {
                    "HA": None,
                    "NA": None,
                }
                by_key[key] = entry
                order.append(key)

            if entry[gene] is not None:
                print(
                    f"WARNING: Multiple distinct {gene} "
                    f"classifications for {key}; keeping first "
                    f"({entry[gene][0]!r} over {seqname!r})",
                    flush=True,
                )
                continue

            entry[gene] = (seqname, row)

    if not order:
        print("No rows for summary; skipping.", flush=True)
        return

    rows_out = []

    for key in order:
        sample, platform, subtype = key
        entry = by_key[key]

        record = {
            "sample": sample,
            "platform": platform,
            "ref_subtype": subtype,
        }

        # Initialize all summary columns to NA.
        for column in SUMMARY_COLUMNS:
            record.setdefault(column, "NA")

        for gene in ("HA", "NA"):
            hit = entry[gene]

            if hit is None:
                continue

            seqname, row = hit

            # Safety net: do not promote a bad-QC hit.
            if is_bad_qc(row.get(QC_COLUMN)):
                print(
                    f"WARNING: Surviving {gene} row for {key} "
                    f"has bad QC ({seqname!r}); writing "
                    f"{gene} columns as NA",
                    flush=True,
                )
                continue

            id_column = (
                "seg4_id" if gene == "HA" else "seg6_id"
            )

            record[id_column] = seqname

            for long_col, summary_col in (
                SUMMARY_GENE_COLUMNS[gene].items()
            ):
                value = row.get(long_col, "NA") or "NA"
                record[summary_col] = value

        rows_out.append(record)

    rows_out.sort(
        key=lambda row: (
            row["sample"],
            row["platform"],
            row["ref_subtype"],
        )
    )

    output_dir = os.path.dirname(output_file)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    with open(output_file, "w", newline="") as out_fh:
        writer = csv.DictWriter(
            out_fh,
            fieldnames=SUMMARY_COLUMNS,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(rows_out)

    print(
        f"Wrote {len(rows_out)} "
        f"(sample, platform, subtype) rows "
        f"-> {output_file}",
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Merge Nextclade HA and NA TSV results into two "
            "long-format tables and one wide summary. "
            "The summary has one row per "
            "(sample, platform, subtype)."
        )
    )

    parser.add_argument(
        "-i",
        "--input-dir",
        required=True,
        help="Directory containing Nextclade TSV result files.",
    )

    parser.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        help=(
            "Prefix for output TSV files. Produces "
            "<prefix>-nextclade_ha.long.tsv, "
            "<prefix>-nextclade_na.long.tsv, and "
            "<prefix>-nextclade_ha_na.short.tsv."
        ),
    )

    args = parser.parse_args()

    pattern = os.path.join(
        args.input_dir,
        "*.nextclade.tsv",
    )

    input_files = sorted(glob.glob(pattern))

    if not input_files:
        raise SystemExit(
            f"No Nextclade TSV result files found in: "
            f"{args.input_dir}"
        )

    print(
        f"Found {len(input_files)} Nextclade TSV files.",
        flush=True,
    )

    # Classify files by gene; keep only HA and NA.
    files_by_gene = {
        "HA": [],
        "NA": [],
    }

    for filename in input_files:
        dataset = get_dataset_name(filename)
        gene = get_gene_from_dataset(dataset)

        if gene not in files_by_gene:
            print(
                f"Skipping (unclassified dataset): {filename}",
                flush=True,
            )
            continue

        files_by_gene[gene].append((filename, gene))

    for gene, files in files_by_gene.items():
        print(
            f"  {gene} files: {len(files)}",
            flush=True,
        )

    wrote_any = False
    long_rows_by_gene = {}

    for gene, files in files_by_gene.items():
        if not files:
            print(
                f"No {gene} input files; skipping output.",
                flush=True,
            )
            continue

        union_columns = collect_union_columns(files)

        print(
            f"Union of {gene} columns: {len(union_columns)}",
            flush=True,
        )

        rows = aggregate_gene(
            files,
            union_columns,
        )

        long_rows_by_gene[gene] = rows

        output_file = LONG_FILENAME_TEMPLATE[gene].format(
            prefix=args.output_prefix,
        )

        write_gene(
            rows,
            union_columns,
            output_file,
        )

        wrote_any = True

    if not wrote_any:
        raise SystemExit(
            "No HA or NA input files were found. Check that the "
            "filenames follow the expected "
            "'all_consensus.nextclade.<dataset>.nextclade.tsv' "
            "pattern."
        )

    ha_rows = long_rows_by_gene.get("HA", [])
    na_rows = long_rows_by_gene.get("NA", [])

    if ha_rows or na_rows:
        summary_file = SUMMARY_FILENAME_TEMPLATE.format(
            prefix=args.output_prefix,
        )

        write_summary(
            ha_rows,
            na_rows,
            summary_file,
        )


if __name__ == "__main__":
    main()
