// subworkflows/local/mapping_nanopore.nf

include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'
include {
    SAMTOOLS_VIEW as SAMTOOLS_VIEW_NANOPORE
} from '../../modules/local/samtools/view/main'

workflow MAPPING_NANOPORE {
    take:
    ch_fastq // channel: [ meta, [reads] ]
    ch_ref_fasta_fai // channel: [ meta, fasta ]

    main:

    ch_fastq
        .join(ch_ref_fasta_fai)
        .multiMap { meta, fastq, fasta, fai ->
            reads: [meta, fastq]
            ref_fasta: [meta, fasta]
        }
        .set { ch_input_minimap2 }

    // 2. Map input variables strictly tracking your process signature
    def bam_format = true
    def bam_index_extension = 'bai'
    // Forces standard .bai index creation
    def cigar_paf_format = false
    def cigar_bam = false

    // 3. Execute the module using your exact input schema
    MINIMAP2_ALIGN(
        ch_input_minimap2.reads,
        ch_input_minimap2.ref_fasta,
        bam_format,
        bam_index_extension,
        cigar_paf_format,
        cigar_bam,
    )

    // 4. Standardize the output tuple to match [meta, bam, bai]
    ch_bam_bai = MINIMAP2_ALIGN.out.bam.join(MINIMAP2_ALIGN.out.index)

    /*
     * Nanopore-specific BAM filtering.
     * Arguments are defined in conf/modules.config.
     */
    ch_input_samtools = ch_bam_bai
        .join(ch_ref_fasta_fai)
        .multiMap { meta, bam, bai, fasta, fai ->
            bam_bai: [meta, bam, bai]
            fasta_fai: [meta, fasta, fai]
        }

    SAMTOOLS_VIEW_NANOPORE(
        ch_input_samtools.bam_bai,
        ch_input_samtools.fasta_fai,
        [[], []],
        [[], []],
        '',
    )

    ch_bam_bai_filtered = SAMTOOLS_VIEW_NANOPORE.out.bam.join(SAMTOOLS_VIEW_NANOPORE.out.csi)

    emit:
    bam_bai = ch_bam_bai_filtered // Emit tracking structural parity with the Illumina mapper
}
