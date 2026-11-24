include {
    BCF_CONSENSUS
} from '../../modules/local/bcftools/bcftools.nf'
include {
    CONSENSUS_REHEADER ;
    CONSENSUS_FILTER_AND_RENAME ;
    CONCAT_CONSENSUS;
    CONSENSUS_REPORT;
    CONCAT_CONSENSUS_REPORT;
} from '../../modules/local/consensus.nf'

include {
    SEQKIT_FX2TAB
} from '../../modules/local/seqkit/fx2tab'



workflow CONSENSUS {
    take:
    fasta_vcf
    per_base_bed
    per_sample_cov

    main:


    //BCFTOOLS_CONSENSUS(vcf_tbi.join(fasta).join(low_depth_bed))
    BCF_CONSENSUS(
        fasta_vcf.join(per_base_bed),
        params.mindepth,
        params.major_allele_fraction
    )

    // CONSENSUS_REHEADER(BCF_CONSENSUS.out.fasta.join(per_sample_cov))
    // CONSENSUS_FILTER(CONSENSUS_REHEADER.out.fasta)
    CONSENSUS_FILTER_AND_RENAME(BCF_CONSENSUS.out.fasta)
    SEQKIT_FX2TAB( CONSENSUS_FILTER_AND_RENAME.out.filtered_fasta, 'txt')
    CONSENSUS_REPORT(SEQKIT_FX2TAB.out.text.join(per_sample_cov).join(CONSENSUS_FILTER_AND_RENAME.out.id_mapping))
    CONSENSUS_REPORT.out.report
        .map { meta, report -> report }
        .collect()
        .set { ch_all_consensus_reports }

    CONCAT_CONSENSUS_REPORT(ch_all_consensus_reports)
    //TODO: rewrite this, working on individual files and then concatenate them together
    // table will contains sample column


    CONCAT_CONSENSUS(
        CONSENSUS_FILTER_AND_RENAME.out.filtered_fasta
            .map { _meta, _fasta -> _fasta }
            .collect()
    )

    all_consensus = CONCAT_CONSENSUS.out.fasta
    report = CONCAT_CONSENSUS_REPORT.out.report

    emit:
    all_consensus
    report
}
