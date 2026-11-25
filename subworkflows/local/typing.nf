include {
    BLAST_BLASTN ;
    PROCESS_BLAST
} from '../../modules/local/blast/blastn'

include {
    CSVTK_ADD_HEADER as CSVTK_ADD_HEADER_BLASTN
} from '../../modules/local/csvtk/add-header/main'


workflow TYPING {
    take:
    fasta // channel: [ meta, [contig] ]

    main:

    BLAST_BLASTN(
        fasta,
        params.typing_db,
    )

    CSVTK_ADD_HEADER_BLASTN(
        BLAST_BLASTN.out.tsv,
        "qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen,qcovs",
    )
    PROCESS_BLAST(CSVTK_ADD_HEADER_BLASTN.out.tsv)
    typing_report = PROCESS_BLAST.out.processed

    emit:
    typing_report
}
