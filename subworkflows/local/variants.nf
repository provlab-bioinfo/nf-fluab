include {
    VARIANTS_ILLUMINA
} from '../../subworkflows/local/variants_illumina'
include {
    VARIANTS_NANOPORE
} from '../../subworkflows/local/variants_nanopore'

include {
    BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_LOWDEPTH
} from '../../modules/nf-core/bedtools/genomecov/main.nf'

include {
    BCF_FILTER;
    BCFTOOLS_STATS;
 } from "../../modules/local/bcftools/bcftools.nf"


workflow VARIANTS {

    take:
        bam_bai //[meta, bam, bai]
        fasta_fai // [meta, fasta, fai]

    main:
        /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        variant calling
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    // 2. DYNAMIC PLATFORM BRANCHING
    bam_bai.branch { meta, _bam, _bai ->
        illumina: meta.platform == 'illumina'
        nanopore: meta.platform == 'nanopore' || meta.platform == 'ont'
    }.set { ch_bam_bai_by_platform }

    VARIANTS_ILLUMINA(ch_bam_bai_by_platform.illumina, fasta_fai)
    VARIANTS_NANOPORE(ch_bam_bai_by_platform.nanopore, fasta_fai)

    // 4. UNIFY ALIGNMENTS (Merge platform streams cleanly)
    ch_unified_vcf = VARIANTS_ILLUMINA.out.vcf.mix(VARIANTS_NANOPORE.out.vcf)

    BCF_FILTER(
        ch_unified_vcf
            .join(fasta_fai)
            .map{ meta, vcf, fasta, _fai ->
                [meta, fasta, vcf]
            },
        params.major_allele_fraction,
        params.minor_allele_fraction

    ) //output [meta, fasta, vcf]

    BCFTOOLS_STATS(BCF_FILTER.out.vcf)

    emit:
        fasta_vcf = BCF_FILTER.out.vcf //[meta, fasta, vcf] filtered vcf

}

