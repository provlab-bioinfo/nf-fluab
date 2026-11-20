process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    // Note: the versions here need to match the versions used in the mulled container below and minimap2/index
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37671219cfd244eb9b33db9345d3543ffd83037419a1c57f4648aace493ec2c2/data' :
        'community.wave.seqera.io/library/minimap2_samtools:b09096fc890429ce' }"
    input:
    tuple val(meta), path(reads)
    path reference
    val sam_format

    output:
    tuple val(meta), path("*.paf"), optional: true, emit: paf
    tuple val(meta), path("*.sam"), optional: true, emit: sam
    tuple val("${task.process}"), val("minimap2"), eval("minimap2 --version"), topic: versions, emit: versions_minimap2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    //def sam_output = sam_format ? "-a | samtools sort | samtools view -@ ${task.cpus} -b -h -o ${prefix}.bam" : "-o ${prefix}.paf"
    def sam_output = sam_format ? "-a -o ${prefix}.sam" : "-o ${prefix}.paf"

    """
    minimap2 \\
        $args \\
        -t $task.cpus \\
        ${reference ?: reads} \\
        $reads \\
        $sam_output


    """
}
