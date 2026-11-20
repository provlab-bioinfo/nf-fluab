include {
    MASH_SCREEN
} from '../../modules/nf-core/mash/screen/main'

include {
    MASHFILTER
} from '../../modules/local/misc'

include {
    SEQKIT_GREP
} from '../../modules/local/seqkit/grep'

include {
    SAMTOOLS_FAIDX
} from '../../modules/local/samtools/faidx'

workflow SEEK_REFERENCES {
    take:
    reads
    msh_db
    fasta_db

    main:

    screen_best = channel.empty()
    ch_msh_db = channel.value(tuple([id: "mshdb"], msh_db))

    // reads.view { "READS: $it" }
    // ch_msh_db.view { "MASH DB: $it" }
    // println "msh_db = ${msh_db}"
    // println "fasta_db = ${fasta_db}"

    MASH_SCREEN(reads, ch_msh_db)

    /*
        If it has hits to multiple types, keep the rows with the best identity for each type.
        For example, if the sample has hits with both H1 and H3, keep the best records of both types.
    */
    MASHFILTER(
        MASH_SCREEN.out.screen.filter { _meta, screen -> screen.countLines() > 0 }
    )
    MASHFILTER.out.screen
        .filter { _meta, screen ->
            screen.countLines() > 1
        }
        .set {
            screen_best
        }

    // [identity, shared-hashes, median-multiplicity, p-value, query-ID, query-comment]
    screen_best
        .map { it ->
            def meta = it[0]
            def arr = it[1].splitCsv(skip: 1, sep: '\t')
            def a = []
            arr.each { n -> a.add(n[4]) }
            [meta, a.join(',')]
        }
        .set { seqids }

    SEQKIT_GREP(seqids, fasta_db)
    ch_ref_fasta = SEQKIT_GREP.out.fasta

    def get_sizes = false

    SAMTOOLS_FAIDX(ch_ref_fasta, get_sizes)
    ch_ref_fasta_fai = ch_ref_fasta.join(SAMTOOLS_FAIDX.out.fai)

    emit:
    screen = screen_best
    ref_fasta = ch_ref_fasta
    ref_fasta_fai = ch_ref_fasta_fai
}
