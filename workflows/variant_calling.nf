nextflow.enable.dsl=2

include { PMDV } from '../modules/variant_calling/pmdv.nf'
include { INDEX_HP    } from '../modules/variant_calling/index_hp.nf'
include { CALL_SV     } from '../modules/variant_calling/call_sv.nf'
include { BIGCLIPPER  } from '../modules/variant_calling/bigclipper.nf'

workflow variant_calling {

    take:
        bam_ch   // tuple(sample, bam, bai)

    main:

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
        bam = bam_ch
}
