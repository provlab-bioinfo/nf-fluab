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
    def model_name = meta.model_name ?: ""

    if (!model_name) {
        error "Unable to determine Clair3 model: meta.model_name is not set."
    }

    // Path to packaged models (inside the container)
    def package_model_base_dir = "/usr/local/bin/models/"

    // Path to user-provided models (staged in the Nextflow work directory)
    def user_model_path = user_model_base_dir ? "\${PWD}/${user_model_base_dir}/${model_name}" : ""

    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Resolve paths using bash inside the container
    PKG_MODEL="${package_model_base_dir}/${model_name}"
    USR_MODEL="${user_model_path}"
    MODEL_PATH=""

    # 1. Check package directory first
    if [ -d "\$PKG_MODEL" ]; then
        MODEL_PATH="\$PKG_MODEL"
    # 2. Check user directory second
    elif [ -n "\$USR_MODEL" ] && [ -d "\$USR_MODEL" ]; then
        MODEL_PATH="\$USR_MODEL"
    # 3. Fallback based on platform
    else
        if [ "${platform}" = "nanopore" ]; then
            MODEL_PATH="${package_model_base_dir}/ont"
        elif [ "${platform}" = "illumina" ]; then
            MODEL_PATH="${package_model_base_dir}/ilmn"
        else
            # If no platform matches, pass original and let Clair3 fail gracefully
            MODEL_PATH="\$PKG_MODEL"
        fi
    fi

    echo "DEBUG: Selected model_name: ${model_name}"
    echo "DEBUG: Resolved model path: \$MODEL_PATH"

    run_clair3.sh \\
        --bam_fn=\${PWD}/${bam} \\
        --ref_fn=\${PWD}/${reference} \\
        --threads=${task.cpus} \\
        --output=. \\
        --platform=${platform} \\
        --model_path=\$MODEL_PATH \\
        ${args}

    # Rename to add prefix
    for file in merge_output.vcf.gz \\
            merge_output.vcf.gz.tbi \\
            phased_merge_output.vcf.gz \\
            phased_merge_output.vcf.gz.tbi \\
            merge_output.gvcf.gz \\
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
