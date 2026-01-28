nextflow.enable.dsl=2

include { DORADO_BASECALLING } from '../modules/basecalling/dorado_basecalling.nf'
include { BASECALLING_REPORT } from '../modules/basecalling/basecalling_report.nf'
include { BAMTOFQ }        from '../modules/alignment/bamtofq.nf'
include { MINIMAP2 }      from '../modules/alignment/minimap2.nf'
include { SAMTOOLS_BAM }  from '../modules/alignment/samtools_bam.nf'
include { SAMTOOLS_TARGET } from '../modules/alignment/samtools_target.nf'
include { FLAGSTAT }      from '../modules/alignment/flagstat.nf'
include { CRAMINO }       from '../modules/alignment/cramino.nf'
include { NANOPLOT }      from '../modules/alignment/nanoplot.nf'
include { COVERAGE }      from '../modules/alignment/coverage.nf'
include { PMDV } from '../modules/variant_calling/pmdv.nf'
include { INDEX_HP    } from '../modules/variant_calling/index_hp.nf'
include { CALL_SV     } from '../modules/variant_calling/call_sv.nf'
include { BIGCLIPPER  } from '../modules/variant_calling/bigclipper.nf'

workflow pipeline {

	take:
   	pod5_ch   // tuple(sample, pod5_dir)

	main:
    	ubam_ch = DORADO_BASECALLING(pod5_ch)
    	report  = BASECALLING_REPORT(ubam_ch)
        
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

        if (!params.skip_PMDV) {

            snv_out = PMDV(bam_ch)

            hp_bam_ch = snv_out[3]  
            phased_snv_vcf = snv_out[2]  
            INDEX_HP(hp_bam_ch)
        } else {
            log.info "Skipping PMDV"
        }

        if (!params.skip_CALL_SV) {
            CALL_SV(bam_ch)
        } else {
            log.info "Skipping CALL_SV"
        }

        BIGCLIPPER(bam_ch)

	emit:
	ubam     = ubam_ch
        bam = bam_ch
}
