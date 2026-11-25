include { FREEBAYES } from '../../modules/nf-core/freebayes/main.nf'
include {
    BCF_FILTER ;
    BCFTOOLS_STATS
} from "../../modules/local/bcftools/bcftools.nf"

workflow VARIANTS_SHORT {
    take:
    bam_bai
    fasta_fai // [meta, bam, bai, fasta, fai]

    main:

    ch_input_freebayes = bam_bai
        .join(fasta_fai)
        .multiMap { meta, bam, bai, fasta, fai ->
            input_1: [meta, bam, bai, [], [], []]
            input_2: [meta, fasta]
            input_3: [meta, fai]
            input_4: [[], []]
            input_5: [[], []]
            input_6: [[], []]
        }

    FREEBAYES(
        ch_input_freebayes.input_1,
        ch_input_freebayes.input_2,
        ch_input_freebayes.input_3,
        ch_input_freebayes.input_4,
        ch_input_freebayes.input_5,
        ch_input_freebayes.input_6,
    )
    BCF_FILTER(
        FREEBAYES.out.vcf.join(fasta_fai).map { meta, vcf, fasta, _fai ->
            [meta, fasta, vcf]
        },
        params.major_allele_fraction,
        params.minor_allele_fraction,
    ) //output [meta, fasta, vcf]

    BCFTOOLS_STATS(BCF_FILTER.out.vcf)

    emit:
    fasta_vcf = BCF_FILTER.out.vcf
}
