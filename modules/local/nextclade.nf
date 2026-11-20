process NEXTCLADE_DATASET_GET {
  tag "$dataset_name"

  conda "bioconda::nextclade=3.23.0"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nextclade:3.23.0--h9ee0642_0' :
        'biocontainers/nextclade:3.23.0--h9ee0642_0' }"

  input:
  val(dataset)

  output:
  tuple val(dataset), path(nextclade_dataset_dir), emit: dir
  path("versions.yml"), emit: versions

  script:
  println "DEBUG dataset.name----------------- = ${dataset.name}"
    println "DEBUG dataset contents:"
    println dataset.text
  dataset_name = dataset.name
  nextclade_dataset_dir = dataset_name.replace('/', '---')
  dataset_tag = dataset.tag ? "--tag ${dataset.tag}" : ""
  """
  nextclade dataset get \\
    --name $dataset_name \\
    $dataset_tag \\
    --output-dir $nextclade_dataset_dir

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/nextclade //')
  END_VERSIONS
  """
}

process NEXTCLADE_RUN {
tag "${meta.id}"

conda "bioconda::nextclade=3.23.0"
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
'https://depot.galaxyproject.org/singularity/nextclade:3.23.0--h9ee0642_0' :
'biocontainers/nextclade:3.23.0--h9ee0642_0' }"

input:
tuple val(meta), path(fasta), val(dataset), path(nextclade_dataset_dir)

output:
tuple val(meta), val(dataset), path(nextclade_tsv), emit: tsv
tuple val(meta), val(dataset), path(nextclade_csv), emit: csv
tuple val(meta), val(dataset), path(nextclade_json), emit: json
tuple val(meta), val(dataset), path(nextclade_warnings), emit: warnings
path("versions.yml"), emit: versions

script:
dataset_name = dataset.name
def prefix = task.ext.prefix ?: "${meta.id}"

nextclade_tsv = "${prefix}.${dataset_name.replace('/', '---')}.nextclade.tsv"
nextclade_json = "${prefix}.${dataset_name.replace('/', '---')}.nextclade.json"
nextclade_csv = "${prefix}.${dataset_name.replace('/', '---')}.nextclade.csv"

nextclade_warnings = "${prefix}.${dataset_name.replace('/', '---')}.nextclade.warnings.log"

"""
nextclade run \\
    --jobs ${task.cpus} \\
    --input-dataset $nextclade_dataset_dir \\
    --output-tsv $nextclade_tsv \\
    --output-json $nextclade_json \\
    --output-csv $nextclade_csv \\
    $fasta \\
    2> $nextclade_warnings

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    nextclade: \$(nextclade --version 2>&1 | sed 's/nextclade //')
END_VERSIONS
"""
}

process AGG_NEXTCLADE_HA_NA {
   tag "aggregate_nextclade"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"


   input:
    path nextclade_csvs

    output:
    path "*-nextclade_ha_na.short.tsv", emit: report
    path "*-nextclade_ha.long.tsv", emit: ha_long
    path "*-nextclade_na.long.tsv", emit: na_long

    script:
    def args = task.ext.args ?: ''

    """
    aggregate_nextclade.py \
        --input-dir . \
        --output-prefix all_consensus
    """
}
