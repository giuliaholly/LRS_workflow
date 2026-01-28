nextflow.enable.dsl=2

include { BAMTOFQ }        from '../modules/alignment/bamtofq.nf'
include { MINIMAP2 }      from '../modules/alignment/minimap2.nf'
include { SAMTOOLS_BAM }  from '../modules/alignment/samtools_bam.nf'
include { SAMTOOLS_TARGET } from '../modules/alignment/samtools_target.nf'
include { FLAGSTAT }      from '../modules/alignment/flagstat.nf'
include { CRAMINO }       from '../modules/alignment/cramino.nf'
include { NANOPLOT }      from '../modules/alignment/nanoplot.nf'
include { COVERAGE }      from '../modules/alignment/coverage.nf'

workflow alignment {

    take:
        ubam_ch

    main:

        fastq_ch = BAMTOFQ(ubam_ch)

        sam_ch      = MINIMAP2(fastq_ch)
        sorted_bam  = SAMTOOLS_BAM(sam_ch)

        bam_ch = params.bed \
            ? SAMTOOLS_TARGET(sorted_bam)
            : sorted_bam

        qc_input = bam_ch.map { sample, bam_file, bam_index ->
            def prefix = params.bed ? "${sample}_target" : sample
            tuple(sample, bam_file, bam_index, prefix)
        }

        if (!params.skip_QC) {
            FLAGSTAT(qc_input)
            CRAMINO(qc_input)
            NANOPLOT(qc_input)
        }

        if (!params.skip_coverage) {
            COVERAGE(qc_input)
        }

    emit:
        bam = bam_ch
}
