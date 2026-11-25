// main.nf or parent workflow


//include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx'
include {
    SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_MAPPING
} from '../../modules/local/samtools/coverage/main'


include {
    BEDTOOLS_GENOMECOV
} from '../../modules/nf-core/bedtools/genomecov/main.nf'

include { MINIBWA_INDEX } from '../../modules/local/minibwa/index'
include { MINIBWA_MAP }   from '../../modules/local/minibwa/map'
include {
    SAMTOOLS_VIEW as SAMTOOLS_VIEW_SHORT;
} from '../../modules/local/samtools/view/main'

workflow MAPPING_SHORT {
    take:
    ch_fastqs       // Queue channel: [ meta, [reads] ]
    ch_ref_fasta_fai    // Queue channel: [ meta, fasta ]

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
    SAMTOOLS_VIEW_SHORT(
        ch_input_samtools.bam_bai,
        ch_input_samtools.fasta_fai,
        [[], []],
        [[], []],
        ''
    )

    ch_bam_bai_filtered = SAMTOOLS_VIEW_SHORT.out.bam.join(SAMTOOLS_VIEW_SHORT.out.csi)


    // 4. UNIFY ALIGNMENTS (Merge platform streams cleanly)
    //ch_unified_bam_bai = MAPPING_ILLUMINA.out.bam_bai.mix(MAPPING_NANOPORE.out.bam_bai)

    // 5. SHARED POST-ALIGNMENT ANALYSIS (Coverage Calculation)
    ch_input_coverage = ch_bam_bai_filtered
        .join(ch_ref_fasta_fai)
        .multiMap { meta, bam, bai, fasta, fai ->
            bam_bai:   [meta, bam, bai]
            fasta_fai: [meta, fasta, fai]
        }

    SAMTOOLS_COVERAGE_MAPPING(
        ch_input_coverage.bam_bai,
        ch_input_coverage.fasta_fai
    )
    def size = []
    def extension = 'bed'
    def sort = 'true'


     //bed file of the low depth regions < params.mindepth
    BEDTOOLS_GENOMECOV(
        ch_bam_bai_filtered.map {
            meta, bam, bai ->
                [meta, bam, 1]
        },
        size,
            extension,
            sort
    )

    emit:
    coverage = SAMTOOLS_COVERAGE_MAPPING.out.coverage
    low_bed =  BEDTOOLS_GENOMECOV.out.genomecov
    bam_bai  =ch_bam_bai_filtered
}
