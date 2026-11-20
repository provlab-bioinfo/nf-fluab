#!/usr/bin/env python3

"""
aggregate_nextclade.py

Merge Nextclade TSV results for influenza HA (segment 4) and
NA (segment 6) into two long-format tables plus one wide
per-(sample, subtype) Nextclade summary table.

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
        segment 4 rows, union header of HA files.

    <prefix>-nextclade_na.long.tsv
        segment 6 rows, union header of NA files.

    <prefix>-nextclade_ha_na.short.tsv
        one row per (sample, subtype), wide HA/NA summary.

    A gene with no input files produces no long-format output file.
    The summary file is written if at least one HA or NA row
    survives filtering.

Long-format output (per gene): one row per surviving input row.

    - The header is the union of all columns across that gene's
      input files, prefixed by source_file and dataset.
    - A file that lacks a given column contributes "NA" for that
      column in its rows.
    - Multiple hits per sample are preserved: a seqName can
      appear in more than one row.
    - No deduplication, no collapsing across files.

Per-(sample, subtype) Nextclade summary output: one row per
(sample, subtype) pair. The output columns are fixed:

    sample
    ref_subtype
    seg4_id
    HA_clade
    HA_proposedSubclade
    HA_subclade
    HA_nextclade_qc.overallStatus
    seg6_id
    NA_clade
    NA_proposedSubclade
    NA_subclade
    NA_nextclade_qc.overallStatus

    sample and ref_subtype are parsed from seqName. Expected
    formats:

        <sample>-seg<4|6>|<subtype>
        <sample>-seg<4|6>|<subtype> ref=<ref_accession>|<ref_subtype>

    Examples:

        flua_i_1-seg4|H3N2
        flua_i_3_co-seg4|H1N1 ref=PV063155|H1N1

    Each summary row carries the HA hit (if any) and the NA hit
    (if any) that match that (sample, subtype). seg4_id / seg6_id
    record the full seqName for that gene on that row, so each
    row points back to the source hits.

    Co-infections (summary only):
        A sample evaluated against more than one subtype (e.g.
        H1N1 and H3N2 reference datasets) produces one summary
        row per (sample, subtype). Each row carries the HA/NA
        hits that matched that subtype.

    Missing values (summary only):
        If a (sample, subtype) has no HA hit, seg4_id and all
        HA_* columns are written as "NA". Same for NA. If a
        field is empty in the input, it is written as "NA".

    Collapsing (summary only):
        Multiple Nextclade results for the same seqName are
        collapsed when they have the same classification. The
        classification key is (clade, proposedSubclade, subclade).
        This handles cases where the same sequence was evaluated
        against multiple reference datasets of the same subtype.
        The first-seen classification wins.

Filtering (applies to both output formats):
    - Rows with a non-empty "errors" field are excluded.
    - Rows with qc.overallStatus == "bad" (case-insensitive,
      whitespace-tolerant) are excluded.
    - As a safety net, the summary also refuses to write a
      "bad" QC value: if one somehow survives, that gene's
      columns on the affected row are written as "NA" and a
      warning is printed.
    - No other QC filtering. qc.overallStatus is carried as a
      column in both outputs.
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


# Gene -> the long-format filename template, relative to the
# user's prefix. "{prefix}" is substituted.
LONG_FILENAME_TEMPLATE = {
    "HA": "{prefix}-nextclade_ha.long.tsv",
    "NA": "{prefix}-nextclade_na.long.tsv",
}


# Filename template for the wide, one-row-per-(sample, subtype)
# summary.
SUMMARY_FILENAME_TEMPLATE = (
    "{prefix}-nextclade_ha_na.short.tsv"
)


PREFIX_COLUMNS = [
    "source_file",
    "dataset",
]


# Columns that define "same classification" for summary
# collapsing.
CLASSIFICATION_COLUMNS = [
    "clade",
    "proposedSubclade",
    "subclade",
]


# The QC column name and the values treated as bad.
QC_COLUMN = "qc.overallStatus"
BAD_QC_VALUE = "bad"


# seqName -> sample / segment / subtype
#
# Accepts both of these forms:
#   flua_i_1-seg4|H3N2
#   flua_i_3_co-seg4|H1N1 ref=PV063155|H1N1
#
# The trailing " ref=...|..." part is optional and ignored.
# The subtype group stops at whitespace or "|", so it captures
# just the query subtype ("H1N1", "H3N2", "Victoria", ...).
SEQNAME_RE = re.compile(
    r"^(?P<sample>.+?)-seg(?P<seg>[0-9]+)\|(?P<subtype>[^|\s]+)"
)


# Per-gene mapping from long-format column -> summary column.
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


# Fixed summary output header. Do not reorder.
SUMMARY_COLUMNS = [
    "sample",
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


def is_bad_qc(value):
    """
    Return True if the value is a "bad" QC status, ignoring case
    and surrounding whitespace.
    """

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
    """
    Determine influenza gene from the Nextclade dataset name.

    Example:
        nextstrain/flu/h3n2/ha/CY163680
        -> HA
    """

    parts = dataset.lower().split("/")

    for part in parts:
        if part in GENE_NAMES:
            return GENE_NAMES[part]

    return None


def read_headers(filename):
    """
    Return the list of fieldnames from a Nextclade TSV file.
    Returns an empty list on read failure.
    """

    try:
        with open(filename, "r", newline="") as in_fh:
            reader = csv.reader(in_fh, delimiter="\t")
            header = next(reader, None)
    except OSError as exc:
        print(f"WARNING: Cannot read {filename}: {exc}", flush=True)
        return []

    return header or []


def collect_union_columns(files):
    """
    Return the union of data columns across the given files,
    preserving first-seen order.

    Columns that collide with PREFIX_COLUMNS are dropped so the
    output header has no duplicates.
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
    """
    Return the tuple that defines "same classification" for a row.
    Used only by the per-(sample, subtype) summary.
    """

    return tuple(
        (row.get(column) or "").strip()
        for column in CLASSIFICATION_COLUMNS
    )


def aggregate_gene(files, union_columns):
    """
    Build one output row per surviving input row for a single gene.

    No collapsing happens here: multiple hits for the same seqName
    are preserved in the long-format tables.

    Returns:
        list of dicts ready for csv.DictWriter.
    """

    rows = []

    total_rows = 0
    kept_rows = 0
    skipped_errors = 0
    skipped_bad_qc = 0

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

            if "errors" not in reader.fieldnames:
                print(
                    f"WARNING: No 'errors' column found: {filename}; "
                    f"error filtering disabled for this file",
                    flush=True,
                )

            if QC_COLUMN not in reader.fieldnames:
                print(
                    f"WARNING: No '{QC_COLUMN}' column found: "
                    f"{filename}; bad-QC filtering disabled for this file",
                    flush=True,
                )

            for row in reader:
                total_rows += 1

                # 1. Drop rows flagged with Nextclade errors.
                errors = (row.get("errors") or "").strip()
                if errors:
                    skipped_errors += 1
                    continue

                # 2. Drop rows with bad overall QC.
                #    Case-insensitive, whitespace-tolerant.
                qc_status = row.get(QC_COLUMN)
                if is_bad_qc(qc_status):
                    skipped_bad_qc += 1
                    continue

                out = {
                    "source_file": os.path.basename(filename),
                    "dataset": dataset,
                }

                for column in union_columns:
                    value = (row.get(column) or "").strip()
                    out[column] = value if value else "NA"

                rows.append(out)
                kept_rows += 1

    print(
        f"  rows read={total_rows}, kept={kept_rows}, "
        f"with errors={skipped_errors}, "
        f"bad QC={skipped_bad_qc}",
        flush=True,
    )

    return rows


def write_gene(rows, union_columns, output_file):
    """
    Write one gene's long-format table:
    PREFIX_COLUMNS + union_columns.

    Every surviving input row is written; nothing is collapsed.
    """

    output_dir = os.path.dirname(output_file)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    fieldnames = PREFIX_COLUMNS + union_columns

    def sort_key(row):
        return (
            row["source_file"],
            row.get("seqName", ""),
        )

    rows_sorted = sorted(rows, key=sort_key)

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
    Build the wide, one-row-per-(sample, subtype) Nextclade
    summary from the long rows.

    ha_rows / na_rows are the lists returned by aggregate_gene(),
    so they already exclude rows with Nextclade errors and rows
    with qc.overallStatus == "bad".

    Safety net: if a row that survived aggregate_gene() still has
    a "bad" QC value in any gene (e.g. because the column was
    missing upstream and the value came from elsewhere, or a
    future Nextclade version emits a different spelling), that
    gene's contribution to the summary row is replaced with
    "NA". The warning is printed so the upstream filter can be
    corrected.

    Co-infections:
        A sample evaluated against more than one subtype produces
        one summary row per (sample, subtype). Each row carries
        the HA/NA hits that matched that subtype; missing genes
        are written as "NA".

    Collapsing:
        Multiple Nextclade results for the same seqName are
        collapsed when they have the same classification
        (clade, proposedSubclade, subclade). The first-seen
        classification wins.
    """

    # per-gene collapse: (seqName, classification) -> first row
    collapsed = {"HA": {}, "NA": {}}
    counts = {"HA": 0, "NA": 0}

    for gene, rows in (("HA", ha_rows), ("NA", na_rows)):
        for row in rows:
            seqname = row.get("seqName", "")
            key = (seqname, classification_key(row))

            if key in collapsed[gene]:
                counts[gene] += 1
                continue

            collapsed[gene][key] = row

        if counts[gene]:
            print(
                f"  summary {gene}: collapsed "
                f"{counts[gene]} duplicate (seqName, classification) "
                f"row(s) from {len(rows)} long-format row(s)",
                flush=True,
            )

    # collect hits per (sample, subtype), per gene
    by_key = {}
    order = []

    for gene in ("HA", "NA"):
        expected_seg = "4" if gene == "HA" else "6"

        for row in collapsed[gene].values():
            seqname = row.get("seqName", "")
            match = SEQNAME_RE.match(seqname)

            if not match:
                print(
                    f"WARNING: cannot parse seqName for summary: "
                    f"{seqname!r}",
                    flush=True,
                )
                continue

            if match.group("seg") != expected_seg:
                print(
                    f"WARNING: segment mismatch for sample "
                    f"{match.group('sample')}: got "
                    f"seg{match.group('seg')} in {gene} file; "
                    f"skipping {seqname!r}",
                    flush=True,
                )
                continue

            sample = match.group("sample")
            subtype = match.group("subtype")
            key = (sample, subtype)

            entry = by_key.get(key)
            if entry is None:
                entry = {"HA": None, "NA": None}
                by_key[key] = entry
                order.append(key)

            if entry[gene] is not None:
                print(
                    f"WARNING: multiple distinct {gene} "
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
        sample, subtype = key
        entry = by_key[key]

        record = {
            "sample": sample,
            "ref_subtype": subtype,
        }

        # Start every summary column at "NA"; fill in below.
        for column in SUMMARY_COLUMNS:
            record.setdefault(column, "NA")

        for gene in ("HA", "NA"):
            hit = entry[gene]
            if hit is None:
                continue

            seqname, row = hit

            # Safety net: never promote a bad-QC hit into the
            # summary. Skip just this gene; keep the other one.
            if is_bad_qc(row.get(QC_COLUMN)):
                print(
                    f"WARNING: surviving {gene} row for {key} "
                    f"has bad QC ({seqname!r}); writing {gene} "
                    f"columns as NA in the summary",
                    flush=True,
                )
                continue

            id_column = "seg4_id" if gene == "HA" else "seg6_id"
            record[id_column] = seqname

            for long_col, summary_col in SUMMARY_GENE_COLUMNS[gene].items():
                record[summary_col] = (
                    row.get(long_col, "NA") or "NA"
                )

        rows_out.append(record)

    rows_out.sort(key=lambda r: (r["sample"], r["ref_subtype"]))

    output_dir = os.path.dirname(output_file)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    with open(output_file, "w", newline="") as out_fh:
        writer = csv.DictWriter(
            out_fh,
            fieldnames=SUMMARY_COLUMNS,
            delimiter="\t",   # change to "," for CSV output
            extrasaction="ignore",
        )

        writer.writeheader()

        for record in rows_out:
            writer.writerow(record)

    print(
        f"Wrote {len(rows_out)} (sample, subtype) rows "
        f"-> {output_file}",
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Merge Nextclade HA and NA TSV results into two long "
            "tables (segment 4 HA, segment 6 NA) plus one wide "
            "per-(sample, subtype) Nextclade report. The long "
            "tables keep every hit; the report has one row per "
            "(sample, subtype) with HA and NA side by side. "
            "Rows with a non-empty 'errors' field or a 'bad' "
            "qc.overallStatus (case-insensitive) are excluded."
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
            "Prefix for the output TSV files. Produces "
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

    print(f"Found {len(input_files)} Nextclade TSV files.", flush=True)

    # Classify files by gene; keep only HA and NA.
    files_by_gene = {"HA": [], "NA": []}

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
            print(f"No {gene} input files; skipping output.", flush=True)
            continue

        union_columns = collect_union_columns(files)

        print(
            f"Union of {gene} columns: {len(union_columns)}",
            flush=True,
        )

        rows = aggregate_gene(files, union_columns)

        long_rows_by_gene[gene] = rows

        output_file = LONG_FILENAME_TEMPLATE[gene].format(
            prefix=args.output_prefix,
        )

        write_gene(rows, union_columns, output_file)

        wrote_any = True

    if not wrote_any:
        raise SystemExit(
            "No HA or NA input files were found. Check that the "
            "filenames follow the expected "
            "'all_consensus.nextclade.<dataset>.nextclade.tsv' pattern."
        )

    ha_rows = long_rows_by_gene.get("HA", [])
    na_rows = long_rows_by_gene.get("NA", [])

    if ha_rows or na_rows:
        summary_file = SUMMARY_FILENAME_TEMPLATE.format(
            prefix=args.output_prefix,
        )
        write_summary(ha_rows, na_rows, summary_file)


if __name__ == "__main__":
    main()
