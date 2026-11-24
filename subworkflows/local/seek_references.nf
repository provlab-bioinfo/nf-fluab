
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

    main:

    // ======================================================
    // Validate reference database
    // ======================================================

    def db_dir = file(params.fluab_reference_db)

    if (!db_dir.exists() || !db_dir.isDirectory()) {
        error "Reference database directory does not exist or is not a directory: ${params.fluab_reference_db}"
    }

    def msh_db = file("${params.fluab_reference_db}/reference.msh")
    def fasta_db = file("${params.fluab_reference_db}/reference.fasta")

    if (!msh_db.exists() || !msh_db.isFile()) {
        error "Mash database file not found: ${msh_db}"
    }

    if (!fasta_db.exists() || !fasta_db.isFile()) {
        error "Reference FASTA file not found: ${fasta_db}"
    }

    log.info "Reference database: ${db_dir.toRealPath()}"
    // ======================================================
    // Mash screening
    // ======================================================

    ch_msh_db = channel.value(tuple([id: "mshdb"], msh_db))

    screen_best = channel.empty()

    MASH_SCREEN(reads, ch_msh_db)

    /*
        If a sample has hits to multiple types, keep the rows
        with the best identity for each type.
    */
    MASHFILTER(
        MASH_SCREEN.out.screen.filter { _meta, screen ->
            screen.countLines() > 0
        }
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
        .map { meta, screen ->
            def arr = screen.splitCsv(skip: 1, sep: '\t')
            def seqids = arr.collect { it[4] }
            [meta, seqids.join(',')]
        }
        .set {
            seqids
        }

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
