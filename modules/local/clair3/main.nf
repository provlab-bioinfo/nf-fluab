process CLAIR3 {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/clair3:2.0.3--py311hbc58adc_0'
        : 'quay.io/biocontainers/clair3:2.0.3--py311hbc58adc_0'}"

    input:
    tuple val(meta), path(bam), path(bai), path(user_model_base_dir), val(platform)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("${prefix}merge_output.vcf.gz"), emit: vcf
    tuple val(meta), path("${prefix}merge_output.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("${prefix}phased_merge_output.vcf.gz"), emit: phased_vcf, optional: true
    tuple val(meta), path("${prefix}phased_merge_output.vcf.gz.tbi"), emit: phased_tbi, optional: true
    tuple val(meta), path("${prefix}merge_output.gvcf.gz"), emit: gvcf, optional: true
    tuple val(meta), path("${prefix}merge_output.gvcf.gz.tbi"), emit: gtbi, optional: true
    tuple val("${task.process}"), val('clair3'), eval('run_clair3.sh --version | sed "s/^Clair3 v//"'), emit: versions_clair3, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def model_name

    if (meta.model_name) {
        model_name = meta.model_name
    }
    else {
        error "Unable to determine Clair3 model: meta.model_name is not set and platform is not nanopore or illumina."
    }
    println "DEBUG: selected model_name = ${model_name}"
    def model_base_dir
    def package_model_base_idr = "/usr/local/bin/models/"

    if (!user_model_base_dir) {
        if (workflow.containerEngine in ['singularity', 'docker', 'podman']) {
            model_base_dir = package_model_base_idr
        }
        else {
            error(
                "Clair3 packaged models are only available in Docker/Singularity/Podman containers. " + "Please use one of these profiles or provide a user_model instead."
            )
        }
    }
    else{
        model_base_dir = user_model_base_dir
    }

    println "DEBUG:  model_base_dir = ${model_base_dir}"

    def model_path ="${package_model_base_idr}/${model_name}"

    if (model_name != 'ont' && model_name != 'ilmn') {

        def model = file("${model_base_dir}/${model_name}")
        println "DEBUG:  model = ${model.toString()}"
        if (model.exists() && model.isDirectory()) {
            model_path = model.toString()
        }
        else if (meta.platform == 'nanopore') {
            model_path = "${package_model_base_idr}/ont"
        }
        else if (meta.platform == 'illumina') {
            model_path = "${package_model_base_idr}/ilmn"
        }
    }


    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    run_clair3.sh \\
        --bam_fn=\${PWD}/${bam} \\
        --ref_fn=\${PWD}/${reference} \\
        --threads=${task.cpus} \\
        --output=. \\
        --platform=${platform} \\
        --model_path=${model_path} \\
        ${args}

    # Rename to add prefix
    for file in merge_output.vcf.gz \
            merge_output.vcf.gz.tbi \
            phased_merge_output.vcf.gz \
            phased_merge_output.vcf.gz.tbi \
            merge_output.gvcf.gz \
            merge_output.gvcf.gz.tbi; do
    if [ -e "\$file" ]; then
        mv "\$file" "${prefix}\$file"
    fi
    done
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}phased_merge_output.vcf.gz
    touch ${prefix}phased_merge_output.vcf.gz.tbi
    echo "" | gzip > ${prefix}merge_output.vcf.gz
    touch ${prefix}merge_output.vcf.gz.tbi
    echo "" | gzip > ${prefix}merge_output.gvcf.gz
    touch ${prefix}merge_output.gvcf.gz.tbi
    """
}
