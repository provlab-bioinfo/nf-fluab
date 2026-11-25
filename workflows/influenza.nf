/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
include { paramsSummaryMap } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_influenza_pipeline'
include {
    SEEK_REFERENCES as SEEK_REFERENCES_SHORT ;
    SEEK_REFERENCES as SEEK_REFERENCES_LONG
} from '../subworkflows/local/seek_references'
include { MAPPING_SHORT } from '../subworkflows/local/mapping_short'
include { MAPPING_LONG } from '../subworkflows/local/mapping_long'
include { VARIANTS_SHORT } from '../subworkflows/local/variants_short'
include { VARIANTS_LONG } from '../subworkflows/local/variants_long'
include {
    CONSENSUS as CONSENSUS_SHORT ;
    CONSENSUS as CONSENSUS_LONG
} from '../subworkflows/local/consensus'
include { TYPING } from '../subworkflows/local/typing'
include { NEXTCLADE } from '../subworkflows/local/nextclade'
include { CREATE_SUMMARY_REPORT } from '../modules/local/misc'
include { VADR_SCAN } from '../modules/local/vadr.nf'

include {
    CONCAT_CONSENSUS ;
    CONCAT_CONSENSUS_REPORT
} from '../modules/local/consensus.nf'

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

    main:
    ch_samplesheet.view { "SAMPLESHEET: ${it}" }

    /*
     * Split by sequencing platform
     */

    ch_illumina = ch_samplesheet.filter { meta, _fastqs ->
        meta.platform == 'illumina'
    }

    ch_nanopore = ch_samplesheet.filter { meta, _fastqs ->
        meta.platform == 'nanopore'
    }

    SEEK_REFERENCES_SHORT(
        ch_illumina
    )
    SEEK_REFERENCES_LONG(
        ch_nanopore
    )

    MAPPING_SHORT(ch_illumina, SEEK_REFERENCES_SHORT.out.ref_fasta_fai)
    MAPPING_LONG(ch_nanopore, SEEK_REFERENCES_LONG.out.ref_fasta_fai)

    VARIANTS_SHORT(MAPPING_SHORT.out.bam_bai, SEEK_REFERENCES_SHORT.out.ref_fasta_fai)
    VARIANTS_LONG(MAPPING_LONG.out.bam_bai, SEEK_REFERENCES_LONG.out.ref_fasta_fai)

    CONSENSUS_SHORT(
        VARIANTS_SHORT.out.fasta_vcf,
        MAPPING_SHORT.out.low_bed,
        MAPPING_SHORT.out.coverage,
    )
    CONSENSUS_LONG(
        VARIANTS_LONG.out.fasta_vcf,
        MAPPING_LONG.out.low_bed,
        MAPPING_LONG.out.coverage,
    )

    /*
     * Combine Illumina and Nanopore consensus FASTAs
     */

    ch_combined_consensus = CONSENSUS_SHORT.out.all_consensus
        .mix(CONSENSUS_LONG.out.all_consensus)
        .collect()

    CONCAT_CONSENSUS(ch_combined_consensus)


    ch_consensus_filtered = CONCAT_CONSENSUS.out.fasta
        .filter { fasta -> fasta.size() > 0 && fasta.countFasta() > 0 }
        .map { fasta ->
            tuple(
                [id: 'all_consensus'],
                fasta,
            )
        }

    VADR_SCAN(ch_consensus_filtered)

    TYPING(ch_consensus_filtered)

    NEXTCLADE(
        ch_consensus_filtered,
        params.nextclade_datasets_csv,
    )


    /*
     * Combine Illumina and Nanopore consensus reports
     */

    ch_combined_consensus_report = CONSENSUS_SHORT.out.report
        .mix(CONSENSUS_LONG.out.report)
        .collect()

    CONCAT_CONSENSUS_REPORT(ch_combined_consensus_report)

    CREATE_SUMMARY_REPORT(
        CONCAT_CONSENSUS.out.fasta,
        CONCAT_CONSENSUS_REPORT.out.report,
        TYPING.out.typing_report.map { _meta, _tsv -> _tsv },
        NEXTCLADE.out.nextclade_report,
        params.flu_db_ver,
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
            storeDir: "${params.outdir}/pipeline_info",
            name: 'influenza_software_' + 'versions.yml',
            sort: true,
            newLine: true,
        )

    emit:
    versions = ch_collated_versions
}
