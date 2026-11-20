process BLAST_BLASTN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::blast=2.17.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.17.0--hb02a186_1':
        'quay.io/biocontainers/blast:2.17.0--hb02a186_1' }"

    input:
    tuple val(meta), path(fasta)
    path  db //fasta reference file

    output:
    tuple val(meta), path('*.typing_output.tsv'), emit: tsv
    tuple val("${task.process}"), val("blastn"), eval("blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//'"), topic: versions, emit: versions_blastn

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    if [[ $fasta == *.gz ]]
    then
        gzip -dc $fasta | blastn -num_threads $task.cpus -subject ${db} $args -out ${prefix}.blastn.tsv
    else
        blastn -num_threads $task.cpus -subject ${db} -query $fasta $args -out ${prefix}.blastn.tsv
    fi

    awk -v OFS='\\t' '!seen[\$1]++' ${prefix}.blastn.tsv > ${prefix}.typing_output.tsv
    """

}
process PROCESS_BLAST {
tag "${meta.id}"

input:
tuple val(meta), path(blast_tsv)

output:
tuple val(meta), path("*.subtype.tsv"), emit: processed

script:
def prefix = task.ext.prefix ?: "${meta.id}"
def output = "${prefix}.subtype.tsv"

"""
awk 'BEGIN {
    OFS="\\t"
    print "qseqid", "subtype", "pident", "qcov", "evalue"
}
NR > 1 {
    split(\$2, a, "~~")
    print \$1, a[3], \$3, \$15, \$11
}' ${blast_tsv} > ${output}
"""
}
