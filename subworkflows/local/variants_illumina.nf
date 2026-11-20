
include { FREEBAYES } from '../../modules/nf-core/freebayes/main.nf'

workflow VARIANTS_ILLUMINA {
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


    emit:
    vcf = FREEBAYES.out.vcf //[meta, vcf.gz]

}
