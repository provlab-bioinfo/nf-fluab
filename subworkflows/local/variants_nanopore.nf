include { CLAIR3 } from '../../modules/local/clair3/main'

workflow VARIANTS_NANOPORE {
    take:
    bam_bai //[meta, bam, bai]
    fasta_fai // [meta, fasta, fai]

    main:

    //bam_bai.view()
    /*
            Clair3 was seeing no read at 2265. Clair3 filters the alignments with the following four flags:
            read unmapped (0x4)
            mate unmapped (0x8)
            not primary alignment (0x100)
            supplementary alignment (0x800)
            Using samtools, the alignments can be filtered with samtools view -F 2316

            Clair3 doesn't call variants in the first 16bp and last 16bp of a sequence because of
            1) algorithmic limit, and
            2) usually degenerated coverage and alignment performance in the head and tail of a sequence
            that makes variant calling unreliable.
            A solution to do the variants calling at the end is to add some 'N' to the tail of your reference genome
            so the end gets out of the tail 16bp limit.
        */

    ch_input = bam_bai
        .join(fasta_fai)
        .multiMap { meta, bam, bai, fasta, fai ->
            input_1: [meta, bam, bai, params.clair3_variant_model, 'ont']
            fasta: [meta, fasta]
            fai: [meta, fai]
        }

    CLAIR3(ch_input.input_1, ch_input.fasta, ch_input.fai)
    vcf = CLAIR3.out.vcf

    emit:
    vcf //[meta, vcf.gz]
}
