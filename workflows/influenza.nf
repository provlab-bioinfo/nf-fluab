/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_influenza_pipeline'



include { READ_QC } from '../subworkflows/local/read_qc'
include { SEEK_REFERENCES } from '../subworkflows/local/seek_references'
include { MAPPING } from '../subworkflows/local/mapping'
include { VARIANTS } from '../subworkflows/local/variants'
include { CONSENSUS } from '../subworkflows/local/consensus'
include { TYPING } from '../subworkflows/local/typing'
include { NEXTCLADE } from '../subworkflows/local/nextclade'
include { CREATE_MASTER_REPORT } from '../modules/local/misc'
include { VADR_SCAN } from '../modules/local/vadr.nf'
/*
references:
    https://www.htslib.org/workflow/wgs-call.html
    https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-00791-w
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow INFLUENZA {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    outdir

    main:
    ch_samplesheet.view { "SAMPLESHEET: $it" }
    READ_QC(ch_samplesheet)
    ch_qcreads = READ_QC.out.qcreads
    //ch_qcreads.view()

    /* SEEK_REFERENCES(
        ch_samplesheet,
        params.flu_db_msh,
        params.flu_db_fasta,
    ) */
    //MAPPING(ch_samplesheet, SEEK_REFERENCES.out.ref_fasta_fai)
    // Find the best matching reference for all samples identically
    SEEK_REFERENCES(
        ch_qcreads,
        params.flu_db_msh,
        params.flu_db_fasta,
    )

    //println "ref = ${SEEK_REFERENCES.out.ref_fasta_fai[1]}"
    MAPPING(ch_qcreads, SEEK_REFERENCES.out.ref_fasta_fai)

    VARIANTS(MAPPING.out.bam_bai, SEEK_REFERENCES.out.ref_fasta_fai)
    CONSENSUS(
        VARIANTS.out.fasta_vcf,
        MAPPING.out.low_bed,
        MAPPING.out.coverage,
    )

    ch_consensus_filtered = CONSENSUS.out.all_consensus
        .filter { fasta -> fasta.size() > 0 && fasta.countFasta() > 0 }
        .map { fasta ->
            tuple(
                [id: "all_consensus"],
                fasta,
            )
        }

    VADR_SCAN(
        ch_consensus_filtered
    )

    TYPING(ch_consensus_filtered)

    NEXTCLADE(
        ch_consensus_filtered,
        params.nextclade_datasets_csv,
    )

    CREATE_MASTER_REPORT(
        CONSENSUS.out.all_consensus,
        CONSENSUS.out.report,
        READ_QC.out.read_qc_report,
        TYPING.out.typing_report.map { _meta, _tsv -> _tsv },
        NEXTCLADE.out.nextclade_report,
        params.flu_db_ver
    )
    // //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(topic_versions.versions_file)
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'influenza_software_' + 'versions.yml',
            sort: true,
            newLine: true,
        )

    emit:
    versions = ch_collated_versions
}
