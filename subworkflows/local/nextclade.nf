nextflow.enable.dsl = 2

include { NEXTCLADE_DATASET_GET ; NEXTCLADE_RUN ; AGG_NEXTCLADE_HA_NA } from '../../modules/local/nextclade.nf'

workflow NEXTCLADE {
    take:
    ch_input_fasta // channel: [sample, path]
    nextclade_datasets_csv // params.nextclade_dataset_csv

    main:
    ch_versions = Channel.empty()

    ch_nextclade_datasets = Channel
        .fromPath(params.nextclade_datasets_csv)
        .splitCsv(header: false, strip: true)
        .map { row ->
            if (row.size() == 0) {
                return null
            }

            def dataset_name = row[0]
            def dataset_tag = row.size() > 1 ? row[1] : null

            return [name: dataset_name, tag: dataset_tag]
        }
        .filter { it != null }

    NEXTCLADE_DATASET_GET(ch_nextclade_datasets)
    ch_versions = ch_versions.mix(NEXTCLADE_DATASET_GET.out.versions)
    ch_datasets = NEXTCLADE_DATASET_GET.out.dir.map { dataset, dataset_dir ->
        // get the actual version tag from the pathogen.json file
        def jsonFile = file("${dataset_dir}/pathogen.json")
        def json = new groovy.json.JsonSlurper().parseText(jsonFile.text)
        dataset.tag = json.version.tag
        return [dataset, dataset_dir]
    }
    NEXTCLADE_RUN(ch_input_fasta.combine(ch_datasets))
    ch_versions = ch_versions.mix(NEXTCLADE_RUN.out.versions)

    AGG_NEXTCLADE_HA_NA(
        NEXTCLADE_RUN.out.tsv.map { _sample, _dataset, _nextclade_tsv -> _nextclade_tsv }.collect()
    )
    nextclade_report = AGG_NEXTCLADE_HA_NA.out.report
    emit:
    nextclade_report
    versions = ch_versions
}
