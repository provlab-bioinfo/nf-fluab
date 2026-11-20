#!/usr/bin/env python

import argparse
from collections import defaultdict
import re
import sys


def readFasta(fasta):
    id2seq = defaultdict(str)
    id2desc = defaultdict(str)

    with open(fasta, "r") as in_f:
        seqid = ""

        for line in in_f:
            if line.startswith(">"):
                full_header = line[1:].strip()

                if full_header:
                    seqid = full_header.split()[0]
                    id2desc[seqid] = full_header
                    id2seq[seqid] = ""
            else:
                id2seq[seqid] += re.sub(r"\s+", "", line)

    return id2seq, id2desc


def filterSeq(id2seq, maxambigs):
    allowed_bases = {"A", "T", "G", "C", "a", "t", "g", "c"}
    excluded_reasons = {}

    for seqid in list(id2seq):
        seq = id2seq[seqid]
        seqlen = len(seq)

        if seqlen == 0:
            excluded_reasons[seqid] = (
                "Reason=Empty_sequence Ambig_Bases=0 "
                "Ambig_Pct=0.00% Total_Length=0"
            )
            del id2seq[seqid]
            continue

        total_atcg_count = len([b for b in seq if b in allowed_bases])
        ambig_count = seqlen - total_atcg_count
        ambig_ratio = ambig_count / seqlen

        if ambig_ratio > maxambigs:
            excluded_reasons[seqid] = (
                f"Reason=Too_many_ambiguities "
                f"Ambig_Bases={ambig_count} "
                f"Ambig_Pct={ambig_ratio:.2%} "
                f"Total_Length={seqlen} "
                f"Limit={maxambigs:.2%}"
            )
            del id2seq[seqid]

    return id2seq, excluded_reasons


def reheader(id2desc, seqPrefix):
    """Build new consensus IDs of the form {prefix}-seg{N}|{SUBTYPE}."""

    id2newid = defaultdict(str)
    id2ref = defaultdict(str)
    id2subtype = defaultdict(str)

    clean_prefix = seqPrefix.replace(".", "_")

    parsed = {}

    for seqid in sorted(id2desc.keys()):
        desc = id2desc[seqid]

        if not desc:
            parsed[seqid] = (None, "", seqid)
            continue

        parts = desc.split("|")

        if len(parts) > 3:
            segid = parts[1].strip()
            subtype = parts[3].strip()
            parsed[seqid] = (segid, subtype, seqid)
        else:
            parsed[seqid] = (None, "", seqid)

    seen = {}

    for seqid, (segid, subtype, ref) in parsed.items():

        if segid is None:
            id2newid[seqid] = f"{clean_prefix}_unmapped_{seqid}"
            id2ref[seqid] = ref
            id2subtype[seqid] = ""
            continue

        tag = subtype if subtype else "unknown"
        newid = f"{clean_prefix}-seg{segid}|{tag}"

        if newid in seen:
            raise ValueError(
                "Duplicate consensus ID would be produced:\n"
                f"  id: {newid}\n"
                f"  sequences: {seen[newid]}, {seqid}"
            )

        seen[newid] = seqid

        id2newid[seqid] = newid
        id2ref[seqid] = ref
        id2subtype[seqid] = subtype

    return id2newid, id2ref, id2subtype


def write_mapping(mapping_out, id2newid):
    """Write new_id -> original_id as a two-column TSV."""

    with open(mapping_out, "w") as out_f:
        out_f.write("new_id\toriginal_id\n")

        for original_id in sorted(id2newid, key=lambda s: id2newid[s]):
            out_f.write(
                f"{id2newid[original_id]}\t{original_id}\n"
            )


def write_fasta(output, sorted_by_new_id, id2ref, id2subtype, id2seq):
    with open(output, "w") as out_f:
        for seqid, newid in sorted_by_new_id:
            ref = id2ref[seqid]
            subtype = id2subtype[seqid]

            if subtype:
                out_f.write(f">{newid} ref={ref}|{subtype}\n")
            else:
                out_f.write(f">{newid} ref={ref}\n")

            out_f.write(f"{id2seq[seqid]}\n")


def write_excluded(excluded_fasta, id2desc, original_id2seq, excluded_reasons):
    with open(excluded_fasta, "w") as ex_f:
        for seqid in sorted(id2desc.keys()):
            if seqid in excluded_reasons:
                stats = excluded_reasons[seqid]
                ex_f.write(f">{id2desc[seqid]} {stats}\n")
                ex_f.write(f"{original_id2seq[seqid]}\n")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Filter consensus sequences by ambiguity ratio, "
            "rename the surviving FASTA headers as "
            "{prefix}-seg{N}|{SUBTYPE}, and write a "
            "new_id -> original_id mapping TSV. Excluded "
            "sequences are written with filter metrics."
        )
    )

    parser.add_argument(
        "-f",
        "--fasta",
        required=True,
        help="FASTA file input",
    )

    parser.add_argument(
        "-p",
        "--prefix",
        required=True,
        help=(
            "Prefix used for the renamed consensus IDs "
            "(e.g. the sample ID) and for the default "
            "excluded output filename"
        ),
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output renamed FASTA file",
    )

    parser.add_argument(
        "-m",
        "--mapping_out",
        help=(
            "Output mapping TSV (new_id -> original_id). "
            "Defaults to {prefix}.id_mapping.tsv"
        ),
    )

    parser.add_argument(
        "--maxambigs",
        type=float,
        default=0.25,
        help="Max ratio of ambiguity bases (default: 0.25)",
    )

    parser.add_argument(
        "-e",
        "--excluded_out",
        help=(
            "Output file path for excluded sequences. "
            "Defaults to {prefix}.excluded.fasta"
        ),
    )

    args = parser.parse_args()

    mapping_out = (
        args.mapping_out
        if args.mapping_out
        else f"{args.prefix}.id_mapping.tsv"
    )

    excluded_fasta = (
        args.excluded_out
        if args.excluded_out
        else f"{args.prefix}.excluded.fasta"
    )

    # 1. Read the original FASTA.
    id2seq, id2desc = readFasta(args.fasta)
    original_id2seq = id2seq.copy()

    # 2. Drop sequences whose ambiguity ratio exceeds the threshold.
    final_id2seq, excluded_reasons = filterSeq(
        id2seq,
        args.maxambigs,
    )

    # 3. Write excluded sequences (original headers + filter metrics).
    write_excluded(
        excluded_fasta,
        id2desc,
        original_id2seq,
        excluded_reasons,
    )

    if excluded_reasons:
        print(
            f"INFO: Wrote {len(excluded_reasons)} excluded "
            f"sequences to '{excluded_fasta}'",
            file=sys.stderr,
        )
    else:
        print(
            f"INFO: No sequences excluded; wrote empty "
            f"'{excluded_fasta}'",
            file=sys.stderr,
        )

    # 4. Rename headers for the surviving sequences only.
    surviving_desc = {
        seqid: id2desc[seqid]
        for seqid in final_id2seq
    }

    id2newid, id2ref, id2subtype = reheader(
        surviving_desc,
        args.prefix,
    )

    # 5. Write the new_id -> original_id mapping TSV.
    write_mapping(
        mapping_out,
        id2newid,
    )

    print(
        f"INFO: Wrote {len(id2newid)} mapping rows "
        f"to '{mapping_out}'",
        file=sys.stderr,
    )

    # 6. Write the renamed FASTA.
    sorted_by_new_id = sorted(
        id2newid.items(),
        key=lambda item: item[1],
    )

    write_fasta(
        args.output,
        sorted_by_new_id,
        id2ref,
        id2subtype,
        final_id2seq,
    )

    print(
        f"INFO: Wrote {len(id2newid)} renamed sequences "
        f"to '{args.output}'",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
