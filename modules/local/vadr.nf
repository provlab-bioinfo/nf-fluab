process VADR_SCAN {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'staphb/vadr:1.7.2' :
        'staphb/vadr:1.7.2' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.vadr.seqstat"),  emit: seqstat
    tuple val(meta), path("${prefix}.*.vadr.pass.tbl"), emit: pass_tbl
    tuple val(meta), path("${prefix}.*.vadr.fail.tbl"), emit: fail_tbl
    tuple val(meta), path("${prefix}.*.vadr.pass.list"), emit: pass_list
    tuple val(meta), path("${prefix}.*.vadr.fail.list"), emit: fail_list
    tuple val(meta), path("${prefix}.*.vadr.pass.fa"),  emit: pass_fa
    tuple val(meta), path("${prefix}.*.vadr.fail.fa"),  emit: fail_fa
    tuple val(meta), path("${prefix}.*.vadr.sqc"),      emit: sqc
    tuple val(meta), path("${prefix}.*.vadr.sqa"),      emit: sqa
    tuple val(meta), path("${prefix}.*.vadr.ftr"),      emit: ftr
    tuple val(meta), path("${prefix}.*.vadr.sgm"),      emit: sgm
    tuple val(meta), path("${prefix}.*.vadr.alt"),      emit: alt
    tuple val(meta), path("${prefix}.*.vadr.dcr"),      emit: dcr
    tuple val(meta), path("${prefix}.*.vadr.rpn"),      emit: rpn
    tuple val(meta), path("${prefix}.*.vadr.mdl"),      emit: mdl
    tuple val(meta), path("${prefix}.*.vadr.alc"),      emit: alc
    tuple val(meta), path("${prefix}.vadr.filelist"), emit: filelist
    tuple val(meta), path("${prefix}.vadr.log"),      emit: log
    tuple val(meta), path("${prefix}.vadr.cmd"),      emit: cmd

    tuple val("${task.process}"), val("vadr"), eval("v-scan.pl -h 2>&1 | sed -n 's/^# VADR \\([0-9.]*\\).*/\\1/p'"), topic: versions, emit: versions_vadr

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    export VADRCONFIGFILE=/opt/vadr/vadr/vadr.config

    v-scan.pl \\
        ${args} \\
        ${fasta} \\
        ${prefix}

    mv ${prefix}/* .
    """
}
