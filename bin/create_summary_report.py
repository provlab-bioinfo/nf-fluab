#!/usr/bin/env python3

"""
create_summary_report.py

Create an influenza analysis summary report for the nf-fluab pipeline.

Inputs:
    1. all_consensus.fasta
    2. report/all_consensus_report.tsv
    3. typing/all_consensus.subtype.tsv
    4. nextclade/ha_na_nextclade_report.tsv

Consensus FASTA headers use:

    >sample-segN|platform|SUBTYPE ref=ACCESSION|SUBTYPE

Example:

    >coinf2i-seg1|illumina|H3N2 ref=PP976202|H3N2

One row represents one sample/platform/reference-subtype combination.

The platform is included in the grouping key so that results from
different sequencing platforms are not combined.

Legacy consensus IDs without a platform field are also supported:

    >sample-segN|SUBTYPE ref=ACCESSION|SUBTYPE

For legacy IDs, platform is reported as "NA".

Reference subtype (ref_subtype)
--------------------------------

ref_subtype is extracted from the ref=ACCESSION|SUBTYPE FASTA header
field. It identifies the reference subtype used for consensus
generation and is used to group consensus segments.

ref_subtype_index is the 1-based ordinal of ref_subtype within each
sample/platform combination, assigned in sorted order.

The per-segment reference accession is not included in the report.

Predicted type or lineage (predicted_type_or_lineage)
------------------------------------------------------

Derived from the per-segment typing results:

    - Influenza A: seg4_typing + seg6_typing
      e.g. H3 + N2 -> H3N2
    - Influenza B: seg4_typing is the lineage
      (Victoria or Yamagata); seg6_typing is ignored.
    - If only one typing result is available, that result is used.
    - If neither is available, the value is "NA".

Per-row genome coverage summary
--------------------------------

avg_genome_completeness and avg_genome_meandepth are arithmetic
means across all eight influenza segments.

Missing or non-numeric segment values count as zero. The denominator
is always eight.

Run-level provenance
--------------------

flu_db_ver and nf_fluab_ver are supplied on the command line and
written into every output row.

Typing results are included for seg4 (HA), seg6 (NA), and seg7 (M).
Nextclade results are included for seg4 (HA) and seg6 (NA).

The Nextclade short report is indexed by seg4_id and seg6_id.
The " ref=..." suffix is stripped before matching consensus IDs.
"""

import argparse
import csv
import re
from collections import defaultdict


SEGMENTS = ("1", "2", "3", "4", "5", "6", "7", "8")

TYPING_SEGMENTS = {"4", "6", "7"}
NEXTCLADE_SEGMENTS = {"4", "6"}

FLU_B_LINEAGES = {"Victoria", "Yamagata"}

CONSENSUS_COLUMNS = [
    "length",
    "numreads",
    "completeness",
    "meandepth",
]

COVERAGE_SUMMARY_METRICS = [
    "completeness",
    "meandepth",
]

TYPING_COLUMNS = [
    "typing",
    "typing_pident",
    "typing_qcov",
    "typing_evalue",
]

NEXTCLADE_GENE_COLUMNS = {
    "HA": {
        "HA_clade": "clade",
        "HA_proposedSubclade": "proposedSubclade",
        "HA_subclade": "subclade",
        "HA_nextclade_qc.overallStatus": "nextclade_qc.overallStatus",
    },
    "NA": {
        "NA_clade": "clade",
        "NA_proposedSubclade": "proposedSubclade",
        "NA_subclade": "subclade",
        "NA_nextclade_qc.overallStatus": "nextclade_qc.overallStatus",
    },
}


def normalize_seq_id(value):
    """Return the bare sequence ID without the reference suffix."""

    value = (value or "").strip()

    if not value or value == "NA":
        return ""

    return value.split(None, 1)[0]


def parse_consensus_id(consensus_id):
    """
    Parse consensus IDs.

    Current format:
        sample-segN|platform|SUBTYPE

    Legacy format:
        sample-segN|SUBTYPE
    """

    match = re.match(
        r"^(?P<sample>.+)-seg(?P<segment>\d+)"
        r"(?P<suffix>.*)$",
        consensus_id,
    )

    if not match:
        raise ValueError(
            f"Unable to parse consensus ID: {consensus_id}"
        )

    sample = match.group("sample")
    segment = match.group("segment")
    suffix = match.group("suffix")

    if not suffix.startswith("|"):
        raise ValueError(
            f"Unable to parse platform/subtype in consensus ID: "
            f"{consensus_id}"
        )

    fields = suffix[1:].split("|")

    if len(fields) == 1:
        platform = "NA"
        subtype = fields[0]

    elif len(fields) == 2:
        platform, subtype = fields

    else:
        raise ValueError(
            f"Unexpected number of fields in consensus ID: "
            f"{consensus_id}"
        )

    if not sample or not subtype:
        raise ValueError(
            f"Missing sample or subtype in consensus ID: "
            f"{consensus_id}"
        )

    return {
        "sample": sample,
        "segment": segment,
        "platform": platform,
        "subtype": subtype,
    }


def read_fasta_headers(fasta):
    """Read and parse consensus FASTA headers."""

    records = {}

    with open(fasta, "r") as in_f:

        for line in in_f:

            if not line.startswith(">"):
                continue

            header = line[1:].strip()
            fields = header.split()

            if not fields:
                continue

            consensus_id = fields[0]

            parsed = parse_consensus_id(consensus_id)

            segment = parsed["segment"]

            if segment not in SEGMENTS:
                raise ValueError(
                    f"Invalid influenza segment '{segment}' "
                    f"in consensus ID: {consensus_id}"
                )

            ref = ""
            reference_subtype = ""

            for field in fields[1:]:

                if field.startswith("ref="):

                    ref_value = field[4:]

                    if "|" in ref_value:
                        ref, reference_subtype = ref_value.split(
                            "|",
                            1,
                        )
                    else:
                        ref = ref_value

                    break

            records[consensus_id] = {
                **parsed,
                "ref": ref,
                "reference_subtype": reference_subtype,
            }

    return records


def read_tsv(path):

    with open(path, "r", newline="") as in_f:

        reader = csv.DictReader(
            in_f,
            delimiter="\t",
        )

        return list(reader)


def read_consensus_report(path):
    """Index consensus QC records by consensus sequence ID."""

    records = {}

    for row in read_tsv(path):

        segment_id = normalize_seq_id(row.get("segment", ""))

        if not segment_id:
            continue

        if segment_id in records:
            raise ValueError(
                f"Multiple consensus report rows found for "
                f"segment: {segment_id}"
            )

        records[segment_id] = {
            column: row.get(column, "NA")
            for column in CONSENSUS_COLUMNS
        }

    return records


def read_typing_report(path):
    """Index typing records by qseqid."""

    records = {}

    for row in read_tsv(path):

        seq_id = normalize_seq_id(row.get("qseqid", ""))

        if seq_id:
            records[seq_id] = row

    return records


def read_nextclade_report(path):
    """
    Read the HA/NA Nextclade short report.

    The report is indexed by seg4_id and seg6_id. The reference
    suffix is removed so the IDs match the consensus FASTA IDs.

    A long-format report with a seqName column is also supported.
    """

    records = {}

    with open(path, "r", newline="") as in_f:

        reader = csv.DictReader(
            in_f,
            delimiter="\t",
        )

        fieldnames = reader.fieldnames or []

        if not any(
            column in fieldnames
            for column in ("seg4_id", "seg6_id", "seqName")
        ):
            raise ValueError(
                f"Nextclade report {path} has none of the expected "
                f"columns: seg4_id, seg6_id, seqName. "
                f"Found: {fieldnames}"
            )

        id_columns = [
            column
            for column in ("seg4_id", "seg6_id")
            if column in fieldnames
        ]

        if not id_columns and "seqName" in fieldnames:
            id_columns = ["seqName"]

        for row_number, row in enumerate(reader, start=2):

            for column in id_columns:

                seq_id = normalize_seq_id(row.get(column, ""))

                if not seq_id:
                    continue

                if seq_id in records:
                    raise ValueError(
                        f"Duplicate nextclade seq ID {seq_id!r} "
                        f"at {path}:{row_number}"
                    )

                records[seq_id] = row

    print(
        f"Read {len(records)} nextclade seq IDs from {path}",
        flush=True,
    )

    return records


def get_consensus_rows(fasta_records):
    """
    Group consensus records by sample, platform, and reference subtype.
    """

    combinations = defaultdict(
        lambda: defaultdict(
            lambda: defaultdict(
                lambda: defaultdict(list)
            )
        )
    )

    for consensus_id, record in fasta_records.items():

        sample = record["sample"]
        platform = record["platform"]
        ref_subtype = record["reference_subtype"]
        segment = record["segment"]

        combinations[sample][platform][ref_subtype][segment].append(
            {
                "consensus_id": consensus_id,
                **record,
            }
        )

    return combinations


def get_single_segment_record(
    sample,
    platform,
    ref_subtype,
    segment,
    records,
):

    if len(records) == 0:
        return None

    if len(records) == 1:
        return records[0]

    consensus_ids = [
        record["consensus_id"]
        for record in records
    ]

    raise ValueError(
        "Multiple consensus sequences found for the same "
        "sample/platform/reference-subtype/segment combination:\n"
        f"  sample: {sample}\n"
        f"  platform: {platform}\n"
        f"  reference subtype: {ref_subtype}\n"
        f"  segment: {segment}\n"
        f"  consensus IDs: {', '.join(consensus_ids)}"
    )


def add_missing_segment(row, segment):

    prefix = f"seg{segment}"

    row[f"{prefix}_id"] = "NA"

    for column in CONSENSUS_COLUMNS:
        row[f"{prefix}_{column}"] = "NA"

    if segment in TYPING_SEGMENTS:

        for column in TYPING_COLUMNS:
            row[f"{prefix}_{column}"] = "NA"

    if segment in NEXTCLADE_SEGMENTS:

        row[f"{prefix}_clade"] = "NA"
        row[f"{prefix}_proposedSubclade"] = "NA"
        row[f"{prefix}_subclade"] = "NA"
        row[f"{prefix}_nextclade_qc.overallStatus"] = "NA"


def add_typing_record(
    row,
    segment,
    consensus_id,
    typing,
):

    prefix = f"seg{segment}"

    typing_record = typing.get(
        consensus_id,
        {},
    )

    row[f"{prefix}_typing"] = typing_record.get(
        "subtype",
        "NA",
    )

    row[f"{prefix}_typing_pident"] = typing_record.get(
        "pident",
        "NA",
    )

    row[f"{prefix}_typing_qcov"] = typing_record.get(
        "qcov",
        "NA",
    )

    row[f"{prefix}_typing_evalue"] = typing_record.get(
        "evalue",
        "NA",
    )


def add_nextclade_record(
    row,
    segment,
    consensus_id,
    nextclade,
):

    prefix = f"seg{segment}"

    nextclade_record = nextclade.get(
        consensus_id,
        {},
    )

    if not nextclade_record:
        print(
            f"WARNING: no nextclade record for {consensus_id!r}",
            flush=True,
        )

    gene = "HA" if segment == "4" else "NA"

    mapping = NEXTCLADE_GENE_COLUMNS[gene]

    for source_column, target_suffix in mapping.items():

        row[f"{prefix}_{target_suffix}"] = (
            nextclade_record.get(source_column, "NA") or "NA"
        )


def add_segment_record(
    row,
    segment,
    record,
    consensus_qc,
    typing,
    nextclade,
):

    prefix = f"seg{segment}"
    consensus_id = record["consensus_id"]

    row[f"{prefix}_id"] = consensus_id

    qc = consensus_qc.get(
        consensus_id,
        {},
    )

    for column in CONSENSUS_COLUMNS:

        row[f"{prefix}_{column}"] = qc.get(
            column,
            "NA",
        )

    if segment in TYPING_SEGMENTS:

        add_typing_record(
            row=row,
            segment=segment,
            consensus_id=consensus_id,
            typing=typing,
        )

    if segment in NEXTCLADE_SEGMENTS:

        add_nextclade_record(
            row=row,
            segment=segment,
            consensus_id=consensus_id,
            nextclade=nextclade,
        )


def derive_predicted_type_or_lineage(row):
    """Derive an Influenza A subtype or Influenza B lineage."""

    ha = (row.get("seg4_typing") or "NA").strip()
    na = (row.get("seg6_typing") or "NA").strip()

    ha_missing = ha in ("", "NA")
    na_missing = na in ("", "NA")

    if not ha_missing and ha in FLU_B_LINEAGES:
        return ha

    if ha_missing and na_missing:
        return "NA"

    if ha_missing:
        return na

    if na_missing:
        return ha

    return f"{ha}{na}"


def add_coverage_summary(row):
    """
    Calculate average completeness and mean depth across all
    eight influenza segments. Missing or non-numeric values count
    as zero.
    """

    segment_count = len(SEGMENTS)

    for metric in COVERAGE_SUMMARY_METRICS:

        total = 0.0

        for segment in SEGMENTS:

            raw = row.get(
                f"seg{segment}_{metric}",
                "NA",
            )

            if raw in (None, "", "NA"):
                continue

            try:
                total += float(raw)
            except (ValueError, TypeError):
                continue

        row[f"avg_genome_{metric}"] = (
            f"{total / segment_count:.2f}"
        )


def build_report(
    fasta_records,
    consensus_qc,
    typing,
    nextclade,
    flu_db_ver,
    nf_fluab_ver,
):

    combinations = get_consensus_rows(fasta_records)

    rows = []

    for sample in sorted(combinations):

        for platform in sorted(combinations[sample]):

            platform_combinations = combinations[sample][platform]

            ref_subtypes = sorted(platform_combinations)

            for ref_subtype_index, ref_subtype in enumerate(
                ref_subtypes,
                start=1,
            ):

                row = {
                    "sample": sample,
                    "platform": platform,
                    "ref_subtype_index": ref_subtype_index,
                    "ref_subtype": ref_subtype or "NA",
                }

                segment_records = platform_combinations[ref_subtype]

                for segment in SEGMENTS:

                    records = segment_records.get(
                        segment,
                        [],
                    )

                    record = get_single_segment_record(
                        sample=sample,
                        platform=platform,
                        ref_subtype=ref_subtype,
                        segment=segment,
                        records=records,
                    )

                    if record is None:

                        add_missing_segment(
                            row,
                            segment,
                        )

                    else:

                        add_segment_record(
                            row=row,
                            segment=segment,
                            record=record,
                            consensus_qc=consensus_qc,
                            typing=typing,
                            nextclade=nextclade,
                        )

                row["predicted_type_or_lineage"] = (
                    derive_predicted_type_or_lineage(row)
                )

                add_coverage_summary(row)

                row["flu_db_ver"] = flu_db_ver
                row["nf_fluab_ver"] = nf_fluab_ver

                rows.append(row)

    return rows


def get_output_columns():

    columns = [
        "sample",
        "platform",
        "ref_subtype_index",
        "ref_subtype",
        "predicted_type_or_lineage",
    ]

    for segment in SEGMENTS:

        prefix = f"seg{segment}"

        columns.extend(
            [
                f"{prefix}_id",
                f"{prefix}_length",
                f"{prefix}_numreads",
                f"{prefix}_completeness",
                f"{prefix}_meandepth",
            ]
        )

        if segment in TYPING_SEGMENTS:

            columns.extend(
                [
                    f"{prefix}_typing",
                    f"{prefix}_typing_pident",
                    f"{prefix}_typing_qcov",
                    f"{prefix}_typing_evalue",
                ]
            )

        if segment in NEXTCLADE_SEGMENTS:

            columns.extend(
                [
                    f"{prefix}_clade",
                    f"{prefix}_proposedSubclade",
                    f"{prefix}_subclade",
                    f"{prefix}_nextclade_qc.overallStatus",
                ]
            )

    columns.extend(
        [
            "avg_genome_completeness",
            "avg_genome_meandepth",
            "flu_db_ver",
            "nf_fluab_ver",
        ]
    )

    return columns


def write_report(rows, output):

    columns = get_output_columns()

    with open(
        output,
        "w",
        newline="",
    ) as out_f:

        writer = csv.DictWriter(
            out_f,
            fieldnames=columns,
            delimiter=",",
            extrasaction="ignore",
        )

        writer.writeheader()

        for row in rows:

            writer.writerow(
                {
                    column: row.get(column, "NA")
                    for column in columns
                }
            )


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Create an influenza analysis summary report "
            "for the nf-fluab pipeline."
        )
    )

    parser.add_argument(
        "-f",
        "--fasta",
        required=True,
        help="All consensus FASTA file.",
    )

    parser.add_argument(
        "--consensus-report",
        required=True,
        help="Consensus QC report TSV.",
    )

    parser.add_argument(
        "--typing-report",
        required=True,
        help="HA/NA/M typing report TSV.",
    )

    parser.add_argument(
        "--nextclade-report",
        required=True,
        help="HA/NA Nextclade short report TSV.",
    )

    parser.add_argument(
        "--flu-db-ver",
        required=True,
        help="Version string of the influenza database used.",
    )

    parser.add_argument(
        "--nf-fluab-ver",
        required=True,
        help="Version string of the nf-fluab pipeline used.",
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output summary report TSV.",
    )

    args = parser.parse_args()

    fasta_records = read_fasta_headers(args.fasta)

    consensus_qc = read_consensus_report(
        args.consensus_report
    )

    typing = read_typing_report(
        args.typing_report
    )

    nextclade = read_nextclade_report(
        args.nextclade_report
    )

    rows = build_report(
        fasta_records=fasta_records,
        consensus_qc=consensus_qc,
        typing=typing,
        nextclade=nextclade,
        flu_db_ver=args.flu_db_ver,
        nf_fluab_ver=args.nf_fluab_ver,
    )

    write_report(
        rows,
        args.output,
    )

    print(
        f"Created summary report: {args.output} "
        f"({len(rows)} rows)",
        flush=True,
    )


if __name__ == "__main__":
    main()
