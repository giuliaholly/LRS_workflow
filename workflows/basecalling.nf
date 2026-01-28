nextflow.enable.dsl=2

include { DORADO_BASECALLING } from '../modules/basecalling/dorado_basecalling.nf'
include { BASECALLING_REPORT } from '../modules/basecalling/basecalling_report.nf'

workflow basecalling {

    take:
    pod5_ch   // tuple(sample, pod5_dir)

    main:
    ubam_ch = DORADO_BASECALLING(pod5_ch)
    report         = BASECALLING_REPORT(ubam_ch)

    emit:
    ubam     = ubam_ch
    report = report
}
