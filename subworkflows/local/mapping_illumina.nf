// subworkflows/local/mapping_illumina.nf

include { MINIBWA_INDEX } from '../../modules/local/minibwa/index'
include { MINIBWA_MAP }   from '../../modules/local/minibwa/map'
include {
    SAMTOOLS_VIEW as SAMTOOLS_VIEW_ILLUMINA;
} from '../../modules/local/samtools/view/main'

workflow MAPPING_ILLUMINA {
    take:
    ch_fastqs       // channel: [ meta, [reads] ]
    ch_ref_fasta_fai    // channel: [ meta, fasta ]

    main:

    ch_fasta = ch_ref_fasta_fai
        .map{
            meta, fasta, _fai ->
                [meta, fasta]
        }

    MINIBWA_INDEX(ch_fasta)

    ch_fastqs
        .join(ch_ref_fasta_fai)
        .join(MINIBWA_INDEX.out.index)
        .multiMap { meta, fastqs, fasta, _fai, index ->
            reads: [meta, fastqs]
            index: [meta, index]
            fasta: [meta, fasta]
        }
        .set { ch_input_minibwa }

    def sort_bam = true // Sent as a proper boolean literal

    MINIBWA_MAP(
        ch_input_minibwa.reads,
        ch_input_minibwa.index,
        ch_input_minibwa.fasta,
        sort_bam
    )

    // Unify the output fields into a paired bam/index channel bundle
    ch_bam_bai = MINIBWA_MAP.out.aligned.join(MINIBWA_MAP.out.index)

    /*
     * Illumina-specific BAM filtering.
     * Arguments are defined in conf/modules.config.
     */
    ch_input_samtools = ch_bam_bai
        .join(ch_ref_fasta_fai)
        .multiMap{
            meta, bam, bai, fasta, fai ->
            bam_bai: [meta, bam, bai]
            fasta_fai: [meta, fasta, fai]
        }
    SAMTOOLS_VIEW_ILLUMINA(
        ch_input_samtools.bam_bai,
        ch_input_samtools.fasta_fai,
        [[], []],
        [[], []],
        ''
    )

    ch_bam_bai_filtered = SAMTOOLS_VIEW_ILLUMINA.out.bam.join(SAMTOOLS_VIEW_ILLUMINA.out.csi)

    emit:
    bam_bai = ch_bam_bai_filtered // emit: [ meta, bam, bai ]
}
