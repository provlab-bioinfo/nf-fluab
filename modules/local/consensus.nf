process CONCAT_CONSENSUS {
    tag "all_samples"

    //publishDir "${params.outdir}/typing", mode: params.publish_dir_mode

    input:
    path(consensus_fastas)

    output:
    path "${prefix}.fasta", emit: fasta

    script:
    prefix = task.ext.prefix ?: "all_consensus"
    """
    cat ${consensus_fastas} > ${prefix}.fasta
    """
}

process CONSENSUS_REHEADER {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(fasta), path(cov)
    output:
    tuple val(meta), path("${prefix}.fa"), emit: fasta
    tuple val(meta), path("*_coverage.tsv"), emit: cov
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed -e "s/Python //g"'), emit: versions_python, topic: versions


    script:

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    consensus_reheader.py \\
        --fasta ${fasta} \\
        --prefix ${meta.id} \\
        --coverage ${cov} \\
        --coverage_out ${prefix}_coverage.tsv \\
        > ${prefix}.fa


    """
}

process CONSENSUS_FILTER_AND_RENAME{
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.fa"), emit: filtered_fasta
    tuple val(meta), path('*excluded.fa'), optional:true, emit: excluded_fasta
    tuple val(meta), path('*id_mapping.tsv'), optional:true, emit: id_mapping
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed -e "s/Python //g"'), emit: versions_python, topic: versions


    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    consensus_filter_and_rename.py \\
        $args \\
        -f ${fasta} \\
        -p ${meta.id} \\
        --platform ${meta.platform} \\
        -e ${prefix}.excluded.fa \\
        -m ${prefix}.id_mapping.tsv \\
        -o ${prefix}.fa

    """
}
process CONSENSUS_REPORT {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(stats), path(coverage), path(id_mapping)

    output:
    tuple val(meta), path("${prefix}.consensus_report.tsv"), emit: report
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed -e "s/Python //g"'), emit: versions_python, topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    consensus_report.py \
        $args \\
        -s ${stats} \\
        -m ${id_mapping} \\
        -c ${coverage} \\
        -o ${prefix}.consensus_report.tsv
    """
}

process CONCAT_CONSENSUS_REPORT {
    tag "all_samples"

    input:
    path(consensus_reports)

    output:
    path "${prefix}.tsv", emit: report

    script:
    prefix = task.ext.prefix ?: "all_consensus_report"

    """
    first=\$(ls ${consensus_reports} | head -1)

    head -n 1 \$first > ${prefix}.tsv

    for report in ${consensus_reports}; do
        tail -n +2 \$report >> ${prefix}.body
    done

    # Sort body by sample (col 1) then segment (col 3).
    # -k1,1 : sample; -k3,3 : segment (string). Use -V for natural seg ordering.
    sort -t \$'\\t' -k1,1 -k3,3V ${prefix}.body >> ${prefix}.tsv

    rm ${prefix}.body
    """
}
