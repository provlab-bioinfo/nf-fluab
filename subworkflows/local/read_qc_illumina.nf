include { FASTP } from '../../modules/nf-core/fastp/main'


workflow READ_QC_ILLUMINA {

    take:
    reads

    main:

    def discard_trimmed_pass = false
    def save_trimmed_fail = false
    def save_merged = false
    def adapter_fasta = params.flu_primers

    FASTP(
        reads.map { _meta, fastqs ->
            [_meta, fastqs, adapter_fasta]
        },
        discard_trimmed_pass,
        save_trimmed_fail,
        save_merged,
    )

    /*
     * Remove zero-size / empty FASTQ files
     */
    trimmed_reads = FASTP.out.reads
        .filter { _meta, _reads ->
            _reads[0].size() > 0 &&
            _reads[0].countFastq() > 0
        }

    emit:
    trimmed_reads
}
