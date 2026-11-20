include { FASTPLONG } from '../../modules/local/fastplong/main'


workflow READ_QC_NANOPORE {

    take:
    reads

    main:

    def discard_trimmed_pass = false
    def save_trimmed_fail = false
    def adapter_fasta = params.flu_primers

    FASTPLONG(
        reads,
        adapter_fasta,
        discard_trimmed_pass,
        save_trimmed_fail,
    )

    /*
     * Remove zero-size / empty FASTQ files
     */
    trimmed_reads = FASTPLONG.out.reads
        .filter { _meta, _reads ->
            _reads.size() > 0 &&
            _reads.countFastq() > 0
        }

    emit:
    trimmed_reads
}
