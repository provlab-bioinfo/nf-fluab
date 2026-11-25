// main.nf or parent workflow


//include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx'
include {
    SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_MAPPING
} from '../../modules/local/samtools/coverage/main'



include{
    MAPPING_ILLUMINA
} from './mapping_illumina.nf'
include{
    MAPPING_NANOPORE
} from './mapping_nanopore.nf'

include {
    BEDTOOLS_GENOMECOV
} from '../../modules/nf-core/bedtools/genomecov/main.nf'


workflow MAPPING {
    take:
    ch_fastqs       // Queue channel: [ meta, [reads] ]
    ch_ref_fasta_fai    // Queue channel: [ meta, fasta ]

    main:

    // 2. DYNAMIC PLATFORM BRANCHING
    ch_fastqs.branch { meta, _reads ->
        illumina: meta.platform == 'illumina'
        nanopore: meta.platform == 'nanopore' || meta.platform == 'ont'
    }.set { ch_reads_by_platform }

    // 3. EXECUTE ISOLATED MAPPING STEP
    MAPPING_ILLUMINA(ch_reads_by_platform.illumina, ch_ref_fasta_fai)
    MAPPING_NANOPORE(ch_reads_by_platform.nanopore, ch_ref_fasta_fai)

    // 4. UNIFY ALIGNMENTS (Merge platform streams cleanly)
    ch_unified_bam_bai = MAPPING_ILLUMINA.out.bam_bai.mix(MAPPING_NANOPORE.out.bam_bai)

    // 5. SHARED POST-ALIGNMENT ANALYSIS (Coverage Calculation)
    ch_input_coverage = ch_unified_bam_bai
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
        ch_unified_bam_bai.map {
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
    bam_bai  = ch_unified_bam_bai
}
