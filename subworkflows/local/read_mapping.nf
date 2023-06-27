//
// Check input samplesheet and get read channels
//

include { BWA_INDEX } from '../../modules/nf-core/bwa/index/main'
include { BWA_MEM } from '../../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_FIXMATE } from '../../modules/nf-core/samtools/fixmate/main'
include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { GATK4_MARKDUPLICATES } from '../../modules/nf-core/gatk4/markduplicates/main'
include { SAMTOOLS_FLAGSTAT } from '../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_COVERAGE } from '../../modules/nf-core/samtools/coverage/main'
include { SAMTOOLS_DEPTH } from '../../modules/nf-core/samtools/depth/main'

workflow READ_MAPPING {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    fasta // channel: [ val(meta), fasta ]

    main:
    ch_versions = Channel.empty()

    // Index reference
    BWA_INDEX( fasta )
    ch_versions = ch_versions.mix( BWA_INDEX.out.versions )

    // Mapping then fixmate and sort bam
    BWA_MEM( reads, BWA_INDEX.out.index, "sort")
        .set { bwa }
    ch_versions = ch_versions.mix( BWA_MEM.out.versions )

    // Markduplication
    GATK4_MARKDUPLICATES( bwa.bam, [], [] )
        .set { bam_markdup }
    ch_versions = ch_versions.mix( GATK4_MARKDUPLICATES.out.versions )

    // Index final bam file
    SAMTOOLS_INDEX( bam_markdup.bam )
        .set { bam_index }
    ch_versions = ch_versions.mix( SAMTOOLS_INDEX.out.versions )

    // Join BAM and BAI channels
    ch_bam = bam_markdup.bam.join( bam_index.bai )

    // Get bam stats: flagstat, coverage, and depth
    SAMTOOLS_DEPTH( bam_markdup.bam )
    SAMTOOLS_FLAGSTAT( ch_bam )
    SAMTOOLS_COVERAGE( ch_bam )

    emit:
    bam = ch_bam                                     // channel: [ val(meta), bam , bai]
    versions = ch_versions                           // channel: [ versions.yml, .. ]
}



