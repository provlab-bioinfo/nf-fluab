

process BCF_CONSENSUS {
  tag "${meta.id}"
  label 'process_medium'

  conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

  input:
  tuple val(meta) , path(fasta), path(vcf), path(per_base_low_bed)
  val(low_coverage)
  val(major_allele_fraction)

  output:
  tuple val(meta), path(consensus), emit: fasta
  tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  consensus    = "${prefix}.bcftools.consensus.fasta"

  """
  # get low coverage depth mask BED file by filtering for regions with less than ${low_coverage}X


  awk '/^>/ {print; next} {gsub(/[RYSWKMBDHVryswkmbdhv]/, "N"); print}' $fasta > ${fasta}.no_ambiguous.fasta

  bcftools filter \\
    -Oz \\
    -o no_low_af_indels.vcf.gz \\
    -e "TYPE != 'snp' && FMT/VAF < ${major_allele_fraction}" \\
    $vcf

  tabix no_low_af_indels.vcf.gz

  bcftools consensus \\
    -f ${fasta}.no_ambiguous.fasta \\
    -m ${per_base_low_bed} \\
    no_low_af_indels.vcf.gz > $consensus


  """
}

process BCF_FILTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta), path(fasta), path(vcf)
    val major_allele_fraction
    val minor_allele_fraction

    output:
    tuple val(meta), path(fasta), path(bcftools_filt_vcf), emit: vcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    bcftools_filt_vcf = "${prefix}.bcftools_filt.vcf"
    """
    bcftools norm \\
        --check-ref w \\
        -Ov \\
        -m- \\
        -f ${fasta} \\
        ${vcf} \\
        > norm.vcf

    bcftools +fill-tags \\
        norm.vcf \\
        -Ov \\
        -o filled.vcf \\
        -- -t all

    bcftools +setGT \\
        filled.vcf \\
        -Ov \\
        -o setGT.major.vcf \\
        -- -t q -n 'c:1/1' -i 'FMT/VAF >= ${major_allele_fraction}'

    bcftools +setGT \\
        setGT.major.vcf \\
        -Ov \\
        -o setGT.minor.vcf \\
        -- -t q  -n 'c:0/1' -i 'FMT/VAF >= ${minor_allele_fraction} && FMT/VAF < ${major_allele_fraction}'

    bcftools +setGT \\
        setGT.minor.vcf \\
        -Ov \\
        -o setGT.final.vcf \\
        -- -t q -n 'c:0/0' -i 'FMT/VAF < ${minor_allele_fraction}'

    bcftools filter \\
        setGT.final.vcf \\
        -e "TYPE != 'snp' && ( (STRLEN(ALT) - STRLEN(REF)) % 3 ) != 0
        || TYPE != 'snp' && ( (STRLEN(ALT) - STRLEN(REF)) % 3 == 0
        && FMT/VAF < ${minor_allele_fraction})" \\
        -Ov \\
        -o ${bcftools_filt_vcf}


    """
}

process BCFTOOLS_STATS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta), path(fasta), path(vcf)

    output:
    path ("*.bcftools_stats.txt"), emit: stats
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    bcftools stats -F ${fasta} ${vcf} > ${prefix}.bcftools_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}

