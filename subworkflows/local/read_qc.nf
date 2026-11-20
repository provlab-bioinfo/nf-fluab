include { READ_QC_ILLUMINA } from './read_qc_illumina'
include { READ_QC_NANOPORE } from './read_qc_nanopore'

include {
    DEACON_FILTER
} from '../../modules/nf-core/deacon/filter/main'

include {
    SEQKIT_STATS as SEQKIT_STATS_INPUT ;
    SEQKIT_STATS as SEQKIT_STATS_TRIMMED ;
    SEQKIT_STATS as SEQKIT_STATS_DEHOSTED
} from '../../modules/local/seqkit/stats/main'

include {
    CSVTK_CONCAT as CSVTK_CONCAT_INPUT_STATS ;
    CSVTK_CONCAT as CSVTK_CONCAT_TRIMMED_STATS ;
    CSVTK_CONCAT as CSVTK_CONCAT_DEHOSTED_STATS
} from '../../modules/nf-core/csvtk/concat/main'

include {
    SEQKIT_STATS_SUMMARY as SEQKIT_STATS_SUMMARY_INPUT;
    SEQKIT_STATS_SUMMARY as SEQKIT_STATS_SUMMARY_TRIMMED;
    SEQKIT_STATS_SUMMARY as SEQKIT_STATS_SUMMARY_DEHOSTED;
    READ_QC_REPORT

 } from "../../modules/local/misc.nf"


workflow READ_QC {
    take:
    ch_samplesheet

    main:

     /*
     * Input read statistics
     */

    SEQKIT_STATS_INPUT(ch_samplesheet)

    CSVTK_CONCAT_INPUT_STATS(
        SEQKIT_STATS_INPUT.out.stats
            .map { _cfg, stats -> stats }
            .collect()
            .map { files ->
                tuple(
                    [id: "seqstats.input"],
                    files,
                )
            },
        "tsv",
        "tsv",
    )
    SEQKIT_STATS_SUMMARY_INPUT(CSVTK_CONCAT_INPUT_STATS.out.csv)
    /*
     * Split by sequencing platform
     */

    ch_illumina = ch_samplesheet.filter { meta, fastqs ->
        meta.platform == 'illumina'
    }

    ch_nanopore = ch_samplesheet.filter { meta, fastqs ->
        meta.platform == 'nanopore'
    }


    /*
     * Platform-specific trimming
     */

    READ_QC_ILLUMINA(ch_illumina)
    READ_QC_NANOPORE(ch_nanopore)


    /*
     * Combine trimmed reads from both platforms
     */

    trimmed_reads = READ_QC_ILLUMINA.out.trimmed_reads.mix(READ_QC_NANOPORE.out.trimmed_reads)
    trimmed_reads.view{}



    /*
    * Trimmed read statistics
    */

    trimmed_reads.view{"trimmed_reads: $it"}
    SEQKIT_STATS_TRIMMED(trimmed_reads)

    CSVTK_CONCAT_TRIMMED_STATS(
        SEQKIT_STATS_TRIMMED.out.stats
            .map { _cfg, stats -> stats }
            .collect().map {
                files ->
                tuple(
                    [id: "seqstats.trimmed"],
                    files,
                )
        }.view(),
        "tsv",
        "tsv",
    )

    SEQKIT_STATS_SUMMARY_TRIMMED(CSVTK_CONCAT_TRIMMED_STATS.out.csv)


    /*
     * Host removal
     */

    trimmed_reads
        .combine(channel.value(params.deacon_refdb))
        .map { meta, fastqs, refdb ->
            [meta, refdb, fastqs]
        }
        .set { ch_input_deacon }

    DEACON_FILTER(ch_input_deacon)

    DEACON_FILTER.out.fastq_filtered
        .view { meta, reads ->
            "DEACON output: meta=${meta}, reads=${reads}"
        }
        .set { ch_deacon_filtered }

    /*
     * Remove zero-size / empty dehosted FASTQ files
     */

    qcreads = DEACON_FILTER.out.fastq_filtered.filter { meta, reads ->
        if (meta.single_end) {
            reads.size() > 0 && reads.countFastq() > 0
        }
        else {
            reads[0].size() > 0 && reads[0].countFastq() > 0
        }
    }


    /*
     * Dehosted read statistics
     */

    SEQKIT_STATS_DEHOSTED(qcreads)

    CSVTK_CONCAT_DEHOSTED_STATS(
        SEQKIT_STATS_DEHOSTED.out.stats
            .map { _cfg, stats -> stats }
            .collect()
            .map { files ->
                tuple(
                    [id: "seqstats.dehosted"],
                    files,
                )
        },
        "tsv",
        "tsv",
    )

    SEQKIT_STATS_SUMMARY_DEHOSTED(CSVTK_CONCAT_DEHOSTED_STATS.out.csv)

    READ_QC_REPORT(
        SEQKIT_STATS_SUMMARY_INPUT.out.summary
            .map{ meta, summary ->
                summary

            },
        SEQKIT_STATS_SUMMARY_TRIMMED.out.summary
            .map{ meta, summary ->
                summary

            },
        SEQKIT_STATS_SUMMARY_DEHOSTED.out.summary
            .map{ meta, summary ->
                summary

            }
    )
    read_qc_report = READ_QC_REPORT.out.report

    emit:
    qcreads
    input_stats = SEQKIT_STATS_INPUT.out.stats
    trimmed_stats = SEQKIT_STATS_TRIMMED.out.stats
    qc_stats = SEQKIT_STATS_DEHOSTED.out.stats
    read_qc_report
}
