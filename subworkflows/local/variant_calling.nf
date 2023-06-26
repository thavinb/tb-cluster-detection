//
// Check input samplesheet and get read channels
//

// original CREATE_SAMPLE_MAP process should be included in GATK4_HAPLOTYPECALLER
include { SAMTOOLS_FAIDX }                  from '../../modules/nf-core/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY }  from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_HAPLOTYPECALLER }           from '../../modules/nf-core/gatk4/haplotypecaller/main'

workflow VARIANT_CALLING {
    take:
    bam_input        // channel: [ val( meta ), bam , bai, interval, drgstrmodel ]
    ch_fasta         // channel: [ val( meta ), fasta ]

    main:
    ch_versions = Channel.empty()

    if (params.dbsnp) { ch_dbsnp = params.dbsnp } else { ch_dbsnp = [] }
    if (params.dbsnp_tbi) { ch_dbsnp_tbi = params.dbsnp_tbi } else { ch_dbsnp_tbi = [] }

    // Index reference
    SAMTOOLS_FAIDX ( ch_fasta )
    ch_versions = ch_versions.mix( SAMTOOLS_FAIDX.out.versions )

    // Create GATK dictionary from reference
    GATK4_CREATESEQUENCEDICTIONARY( ch_fasta )
    ch_versions = ch_versions.mix( GATK4_CREATESEQUENCEDICTIONARY.out.versions )

    // Variant Calling
    // TODO: make sample map from vcf output.
    GATK4_HAPLOTYPECALLER (
        bam_input,
        ch_fasta[1],
        SAMTOOLS_FAIDX.out.fai
            .map { it[1] },
        GATK4_CREATESEQUENCEDICTIONARY.out.dict
            .map { it[1] },
        ch_dbsnp,
        ch_dbsnp_tbi
    )
    ch_versions = ch_versions.mix( GATK4_HAPLOTYPECALLER.out.versions )

    emit:
    vcf = GATK4_HAPLOTYPECALLER.out.vcf                              // channel: [ val(meta), vcf ]
    tbi = GATK4_HAPLOTYPECALLER.out.tbi                              // channel: [ val(meta), tbi ]
    versions = ch_versions                                           // channel: [ versions.yml ]
}

