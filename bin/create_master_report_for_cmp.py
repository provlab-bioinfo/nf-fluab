#!/usr/bin/env python3

"""
create_master_report.py

Create a master influenza analysis report by combining:

    1. all_consensus.fasta
    2. report/all_consensus_report.tsv
    3. report/read_qc_report.tsv
    4. typing/all_consensus.subtype.tsv
    5. nextclade/ha_na_nextclade_report.tsv

Design
------

Consensus FASTA headers use:

    >sample-segN|SUBTYPE ref=ACCESSION|SUBTYPE

Examples:

    >coinf2i-seg1|H3N2 ref=PP976202|H3N2
    >coinf2i-seg4|H1N1 ref=PV063155|H1N1

The |SUBTYPE suffix identifies the reference subtype used for that
consensus sequence. Within a (sample, segment), each subtype appears
at most once.

Segment numbers follow the influenza convention (seg1..seg8, in
decreasing segment length). The protein identity of a given segment
number depends on the influenza type: for Influenza A, seg1=PB2 and
seg2=PB1; for Influenza B, seg1=PB1 and seg2=PB2. The other segments
(PA, HA, NP, NA, M, NS) are the same in both. This script does not
attach gene names to segment numbers; output columns are keyed by
segment number only.

One row represents one sample/reference-subtype combination.

Reference subtype (ref_subtype)
-------------------------------

ref_subtype is the |SUBTYPE from the FASTA header (i.e. the
reference subtype used to assemble the consensus). It is a
grouping signal based on the seekref reference and should not be
interpreted as definitive biological phasing.

ref_subtype_index is the 1-based ordinal of ref_subtype within the
sample, assigned in sorted order of the ref_subtype strings. It is
a stable identifier for the (sample, ref_subtype) row.

The per-segment reference accession is not carried into the
master report; the segN_ref column of earlier revisions has been
removed. The accession remains available in the FASTA headers and
in the long nextclade tables if needed.

Predicted type or lineage (predicted_type_or_lineage)
-----------------------------------------------------

predicted_type_or_lineage is derived from the per-segment typing
results and holds either an Influenza A subtype or an Influenza B
lineage:

    - Influenza A: subtype = seg4_typing (HA) +
      seg6_typing (NA), e.g. "H3" + "N2" -> "H3N2".
      If only one half is present, that half is used as-is
      (e.g. "H3" or "N2").

    - Influenza B: lineage = seg4_typing (Victoria or
      Yamagata). Influenza B has no HA/NA subtype split, so the
      HA typing *is* the lineage and seg6 typing is not
      consulted.

    - If neither is available, predicted_type_or_lineage = "NA".

The name "predicted_type_or_lineage" is used because the column
holds different concepts for the two influenza types (a subtype
for Flu A, a lineage for Flu B) and because it is predicted from
the typing results rather than taken from the reference.

Per-row genome coverage summary
(avg_genome_completeness, avg_genome_meandepth)
----------------------------------------------

avg_genome_completeness and avg_genome_meandepth are the
arithmetic means of segN_completeness and segN_meandepth over
all influenza segments. The denominator is always the total
number of segments (8), so missing segments count as zero and
dilute the average. A non-numeric per-segment value is treated
as zero.

These two columns are appended at the very end of the report,
after all per-segment columns.

Run-level provenance (flu_db_ver, nf_fluab_ver)
-----------------------------------------------

flu_db_ver and nf_fluab_ver are run-level constants supplied on
the command line:

    --flu-db-ver    version of the influenza database used
    --nf-fluab-ver  version of the nf-fluab pipeline used

Both values are written into every output row so that any single
row of the report is self-contained for provenance. They are the
last two columns of the report.

Typing results are included for seg4 (HA), seg6 (NA), seg7 (M).
Nextclade results are included for seg4 (HA) and seg6 (NA).

The nextclade report pairs HA and NA in the same row by reference
subtype. read_nextclade_report() indexes those rows by both
seg4_id and seg6_id so that each consensus ID can look up its own
gene's clade columns.

The seg4_id / seg6_id values in the nextclade report carry a
" ref=ACCESSION|SUBTYPE" suffix (e.g.
"flua_i_1-seg4|H3N2 ref=PV299295|H3N2"). The consensus IDs used
elsewhere in this script (FASTA headers, consensus QC report,
typing report) are the bare form without that suffix
("flua_i_1-seg4|H3N2"). read_nextclade_report() strips the suffix
by taking the first whitespace-delimited token, so the two forms
match.

The short nextclade report carries gene-prefixed classification
columns (HA_clade, HA_proposedSubclade, HA_subclade,
HA_nextclade_qc.overallStatus, and the NA equivalents). The
NEXTCLADE_GENE_COLUMNS mapping translates those prefixed source
columns into the unprefixed suffixes used by the master report's
segN_* columns.
"""

import argparse
import csv
import re
from collections import defaultdict


# Influenza segment numbers, in output column order.
# The protein identity of a segment depends on the influenza
# type (Flu A vs Flu B), so no gene name is attached here.
SEGMENTS = ("1", "2", "3", "4", "5", "6", "7", "8")


TYPING_SEGMENTS = {"4", "6", "7"}
NEXTCLADE_SEGMENTS = {"4", "6"}


# Influenza B lineages. Flu B has no HA/NA subtype split; the
# lineage is reported by the HA (seg4) typing.
FLU_B_LINEAGES = {"Victoria", "Yamagata"}


CONSENSUS_COLUMNS = [
    "length",
    "numreads",
    "completeness",
    "meandepth",
]


# Metrics averaged per row over all 8 segments.
COVERAGE_SUMMARY_METRICS = [
    "completeness",
    "meandepth",
]


READ_QC_COLUMNS = [
    "read_type",
    "input_num_reads",
    "input_total_bases",
    "input_avg_len",
    "trimmed_num_reads",
    "trimmed_total_bases",
    "trimmed_avg_len",
    "trimmed_pct",
    "dehosted_num_reads",
    "dehosted_total_bases",
    "dehosted_avg_len",
    "dehosted_pct_of_trimmed",
    "dehosted_pct_of_input",
]


TYPING_COLUMNS = [
    "typing",
    "typing_pident",
    "typing_qcov",
    "typing_evalue",
]


# Per-gene mapping from nextclade short-report column ->
# master-report suffix.
#
# The short report produced by aggregate_nextclade.py prefixes
# every classification column with the gene name:
#
#     HA_clade, HA_proposedSubclade, HA_subclade,
#     HA_nextclade_qc.overallStatus,
#     NA_clade, NA_proposedSubclade, NA_subclade,
#     NA_nextclade_qc.overallStatus
#
# The master report writes them as seg4_<suffix> / seg6_<suffix>,
# so the mapping below strips the gene prefix by mapping the
# prefixed source name to the unprefixed target suffix.
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
    """
    Return the bare seqName without the " ref=...|..." suffix.

    The seqName is the first whitespace-delimited token, so this
    handles any whitespace or invisible characters before the
    suffix, and both "ref=ACCESSION|SUBTYPE" and "ref=ACCESSION".

    Returns "" for empty or "NA" values.
    """

    value = (value or "").strip()

    if not value or value == "NA":
        return ""

    return value.split(None, 1)[0]


def read_fasta_headers(fasta):
    """
    Read consensus FASTA headers.

    Expected format:

        >sample-segN|SUBTYPE ref=ACCESSION|SUBTYPE

    The |SUBTYPE suffix is optional (legacy untagged IDs still parse).
    """

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

            match = re.match(
                r"^(?P<sample>.+)-seg(?P<segment>\d+)"
                r"(?:\|[^|]+)?$",
                consensus_id,
            )

            if not match:
                raise ValueError(
                    f"Unable to parse consensus ID: {consensus_id}"
                )

            sample = match.group("sample")
            segment = match.group("segment")

            if segment not in SEGMENTS:

                raise ValueError(
                    f"Invalid influenza segment '{segment}' "
                    f"in consensus ID: {consensus_id}"
                )

            records[consensus_id] = {
                "sample": sample,
                "segment": segment,
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
    """
    Read the consensus report TSV.

    The report column that holds the consensus sequence name is
    called "segment" and contains the full consensus ID
    (e.g. flua_i_4_co-seg1|H3N2). We key on that value so that
    FASTA header consensus IDs can be matched to their QC row.
    """

    records = {}

    for row in read_tsv(path):

        segment_id = row["segment"]

        if segment_id in records:
            raise ValueError(
                f"Multiple consensus report rows found for segment: {segment_id}"
            )

        records[segment_id] = {
            column: row.get(column, "NA")
            for column in CONSENSUS_COLUMNS
        }

    return records


def read_read_qc_report(path):
    records = {}

    for row in read_tsv(path):

        sample = row["sample"]

        if sample in records:
            raise ValueError(
                f"Multiple read QC records found for sample: {sample}"
            )

        records[sample] = row

    return records


def read_typing_report(path):
    records = {}

    for row in read_tsv(path):
        records[row["qseqid"]] = row

    return records


def read_nextclade_report(path):
    """
    Read the HA/NA Nextclade short report.

    Expected input: the wide per-(sample, reference subtype) report
    produced by aggregate_nextclade.py, with columns

        seg4_id, HA_clade, HA_proposedSubclade, HA_subclade,
        HA_nextclade_qc.overallStatus,
        seg6_id, NA_clade, NA_proposedSubclade, NA_subclade,
        NA_nextclade_qc.overallStatus

    seg4_id / seg6_id carry a " ref=ACCESSION|SUBTYPE" suffix; the
    suffix is stripped so lookups by bare consensus ID succeed.

    Fallback: if the input has a "seqName" column instead of
    seg4_id / seg6_id (e.g. one of the long-format tables), each
    row is indexed under its own seqName.

    Duplicate bare IDs are a hard error.
    """

    records = {}

    with open(path, "r", newline="") as in_f:

        reader = csv.DictReader(in_f, delimiter="\t")

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

        for row_number, row in enumerate(reader, start=2):

            id_columns = [
                column
                for column in ("seg4_id", "seg6_id")
                if column in fieldnames
            ]

            if not id_columns and "seqName" in fieldnames:
                id_columns = ["seqName"]

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

    combinations = defaultdict(
        lambda: defaultdict(
            lambda: defaultdict(list)
        )
    )

    for consensus_id, record in fasta_records.items():

        sample = record["sample"]
        ref_subtype = record["reference_subtype"]
        segment = record["segment"]

        combinations[sample][ref_subtype][segment].append(
            {
                "consensus_id": consensus_id,
                **record,
            }
        )

    return combinations


def get_single_segment_record(
    sample,
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
        "sample/reference-subtype/segment combination:\n"
        f"  sample: {sample}\n"
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

    mapping = NEXTCLADE_GENE_COLUMNS["HA" if segment == "4" else "NA"]

    for source_column, target_suffix in mapping.items():

        row[f"{prefix}_{target_suffix}"] = nextclade_record.get(
            source_column,
            "NA",
        ) or "NA"


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
    """
    Return the typing-derived predicted type or lineage for a row.

    The value is either an Influenza A subtype ("H3N2", "H1N1",
    ...) or an Influenza B lineage ("Victoria", "Yamagata").

    Rules:

        1. If seg4_typing is a Flu B lineage
           (Victoria / Yamagata), that value is the type;
           seg6_typing is ignored.
        2. Otherwise, if both seg4 and seg6 typing are present,
           concatenate them: "H3" + "N2" -> "H3N2".
        3. If only seg4 typing is present, use it: "H3".
        4. If only seg6 typing is present, use it: "N1".
        5. If neither is present, return "NA".
    """

    ha = (row.get("seg4_typing") or "NA").strip()
    na = (row.get("seg6_typing") or "NA").strip()

    ha_missing = ha in ("", "NA")
    na_missing = na in ("", "NA")

    # Flu B: seg4 typing is the lineage; ignore seg6.
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
    Add avg_genome_completeness and avg_genome_meandepth to the row.

    Each average is the sum of the per-segment values divided by
    the total number of influenza segments (len(SEGMENTS), i.e. 8).
    Missing segments count as zero, so a row with only some of the
    8 segments present will have its averages diluted accordingly.

    Non-numeric per-segment values are treated as zero.
    """

    segment_count = len(SEGMENTS)

    for metric in COVERAGE_SUMMARY_METRICS:

        total = 0.0

        for segment in SEGMENTS:

            raw = row.get(f"seg{segment}_{metric}", "NA")

            if raw in (None, "", "NA"):
                continue

            try:
                total += float(raw)
            except ValueError:
                continue

        row[f"avg_genome_{metric}"] = (
            f"{total / segment_count:.2f}"
        )


def build_report(
    fasta_records,
    consensus_qc,
    read_qc,
    typing,
    nextclade,
    flu_db_ver,
    nf_fluab_ver,
):

    combinations = get_consensus_rows(
        fasta_records
    )

    rows = []

    # read_qc is the master sample list. This ensures that samples
    # without any consensus sequence are still included.
    for sample in sorted(read_qc):

        sample_qc = read_qc[sample]

        sample_combinations = combinations.get(
            sample,
            {},
        )

        # Keep samples that have no consensus.
        if not sample_combinations:

            row = {
                "sample": sample,
                "ref_subtype_index": 1,
                "ref_subtype": "NA",
                "predicted_type_or_lineage": "NA",
            }

            for segment in SEGMENTS:

                add_missing_segment(
                    row,
                    segment,
                )

            add_coverage_summary(row)

            row["flu_db_ver"] = flu_db_ver
            row["nf_fluab_ver"] = nf_fluab_ver

            for column in READ_QC_COLUMNS:

                row[column] = sample_qc.get(
                    column,
                    "NA",
                )

            rows.append(row)
            continue

        ref_subtypes = sorted(
            sample_combinations
        )

        for ref_subtype_index, ref_subtype in enumerate(
            ref_subtypes,
            start=1,
        ):

            row = {
                "sample": sample,
                "ref_subtype_index": ref_subtype_index,
                "ref_subtype": ref_subtype or "NA",
            }

            segment_records = sample_combinations[
                ref_subtype
            ]

            for segment in SEGMENTS:

                records = segment_records.get(
                    segment,
                    [],
                )

                record = get_single_segment_record(
                    sample=sample,
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

            # Derive `predicted_type_or_lineage` and coverage
            # summary only after all segments (and thus all typing
            # and QC fields) have been populated.
            row["predicted_type_or_lineage"] = (
                derive_predicted_type_or_lineage(row)
            )
            add_coverage_summary(row)

            row["flu_db_ver"] = flu_db_ver
            row["nf_fluab_ver"] = nf_fluab_ver

            for column in READ_QC_COLUMNS:

                row[column] = sample_qc.get(
                    column,
                    "NA",
                )

            rows.append(row)

    return rows


def get_output_columns():

    columns = [
        "sample",
    ]

    columns.extend(
        READ_QC_COLUMNS
    )

    columns.extend(
        [
            "ref_subtype_index",
            "ref_subtype",
            "predicted_type_or_lineage",
        ]
    )

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

    # Row-level genome coverage summaries, then run-level
    # provenance. Appended at the very end so they read as a
    # roll-up of the per-segment columns above.
    columns.extend(
        [
            "avg_genome_completeness",
            "avg_genome_meandepth",
            "flu_db_ver",
            "nf_fluab_ver",
        ]
    )

    return columns


def write_report(
    rows,
    output,
):

    columns = get_output_columns()

    with open(
        output,
        "w",
        newline="",
    ) as out_f:

        writer = csv.DictWriter(
            out_f,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writeheader()

        for row in rows:

            writer.writerow(
                {
                    column: row.get(
                        column,
                        "NA",
                    )
                    for column in columns
                }
            )


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Create a master influenza analysis report using read QC "
            "as the sample-level master."
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
        "--read-qc-report",
        required=True,
        help="Read QC report TSV.",
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
        help=(
            "Version string of the influenza database used. "
            "Written to every row as flu_db_ver."
        ),
    )

    parser.add_argument(
        "--nf-fluab-ver",
        required=True,
        help=(
            "Version string of the nf-fluab pipeline used. "
            "Written to every row as nf_fluab_ver."
        ),
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output master report TSV.",
    )

    args = parser.parse_args()

    fasta_records = read_fasta_headers(
        args.fasta
    )

    consensus_qc = read_consensus_report(
        args.consensus_report
    )

    read_qc = read_read_qc_report(
        args.read_qc_report
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
        read_qc=read_qc,
        typing=typing,
        nextclade=nextclade,
        flu_db_ver=args.flu_db_ver,
        nf_fluab_ver=args.nf_fluab_ver,
    )

    write_report(
        rows,
        args.output,
    )


if __name__ == "__main__":
    main()
