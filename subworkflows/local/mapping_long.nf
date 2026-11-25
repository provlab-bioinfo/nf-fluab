// main.nf or parent workflow


//include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx'
include {
    SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_MAPPING
} from '../../modules/local/samtools/coverage/main'

include {
    BEDTOOLS_GENOMECOV
} from '../../modules/nf-core/bedtools/genomecov/main.nf'

include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'
include {
    SAMTOOLS_VIEW as SAMTOOLS_VIEW_LONG
} from '../../modules/local/samtools/view/main'

workflow MAPPING_LONG {
    take:
    ch_fastqs       // Queue channel: [ meta, [reads] ]
    ch_ref_fasta_fai    // Queue channel: [ meta, fasta ]

    main:

    // 3. EXECUTE ISOLATED MAPPING STEP

    //MAPPING_NANOPORE(ch_reads_by_platform.nanopore, ch_ref_fasta_fai)
    ch_fastqs
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

    SAMTOOLS_VIEW_LONG(
        ch_input_samtools.bam_bai,
        ch_input_samtools.fasta_fai,
        [[], []],
        [[], []],
        '',
    )

    ch_bam_bai_filtered = SAMTOOLS_VIEW_LONG.out.bam.join(SAMTOOLS_VIEW_LONG.out.csi)

    // 4. UNIFY ALIGNMENTS (Merge platform streams cleanly)
    //ch_unified_bam_bai = MAPPING_ILLUMINA.out.bam_bai.mix(MAPPING_NANOPORE.out.bam_bai)

    // 5. SHARED POST-ALIGNMENT ANALYSIS (Coverage Calculation)
    ch_input_coverage =  ch_bam_bai_filtered
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
    bam_bai  =  ch_bam_bai_filtered
}
