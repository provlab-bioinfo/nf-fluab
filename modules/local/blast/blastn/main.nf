process BLAST_BLASTN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0c/0c86cbb145786bf5c24ea7fb13448da5f7d5cd124fd4403c1da5bc8fc60c2588/data':
        'community.wave.seqera.io/library/blast:2.17.0--d4fb881691596759' }"

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
        print "sample", "platform", "seqid", "subtype", "pident", "qcov", "evalue"
    }
    NR > 1 {
        split(\$1, q, "|")
        split(q[1], s, "-seg")

        if (length(q) == 3 && length(s) == 2) {
            sample = s[1]
            platform = q[2]
            seqid = \$1
            subtype = q[3]

            print sample, platform, seqid, subtype, \$3, \$15, \$11
        }
    }' ${blast_tsv} > ${output}
    """
}

