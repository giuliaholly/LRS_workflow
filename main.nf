nextflow.enable.dsl = 2

if (params.help) {
    log.info "Use --samplesheet samplesheet.tsv"
    exit 0
}

if (!params.samplesheet) {
    error "You must provide --samplesheet"
}

/*
 * Import modules
 */
include { BASECALLING; BASECALLING_REPORT } from './modules/basecalling.nf'
include { BAMTOFQ; MINIMAP2; SAMTOOLS_BAM; SAMTOOLS_TARGET; FLAGSTAT; CRAMINO; NANOPLOT; COVERAGE } from './modules/alignment.nf'
include { SNVs_PEPPER; INDEX_HP; CALL_SV; BIGCLIPPER } from './modules/variant_calling.nf'

workflow {

    /*
     * Parse samplesheet
     */
samples_ch = Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true)  // usa default = ','
    .map { row ->
        def sample = row.sample
        if (!sample) error "Missing sample name in samplesheet row: ${row}"

        if (!params.skip_basecalling && !params.skip_alignment) {
            if (!row.pod5 || row.pod5 == '.') error "Missing pod5 for sample ${sample}"
            return tuple(sample, file(row.pod5), 'POD5')
        } else if (params.skip_basecalling && !params.skip_alignment) {
            if (!row.fastq || row.fastq == '.') error "Missing FASTQ for sample ${sample}"
            return tuple(sample, file(row.fastq), 'FASTQ')
        } else if (params.skip_alignment) {
            if (!row.bam || row.bam == '.') error "Missing BAM for sample ${sample}"
            return tuple(sample, file(row.bam), file(row.bam+'.bai'), 'BAM')
        } else {
            error "Invalid samplesheet configuration for sample ${sample}"
        }
    }

    /*
     * STAGE 1
     */
if (!params.skip_basecalling && !params.skip_alignment) {

    pod5_ch = samples_ch
        .filter { it[2] == 'POD5' }
        .map { sample, pod5, _ -> tuple(sample, pod5) }

    basecalled_ch = BASECALLING(pod5_ch)
    BASECALLING_REPORT(basecalled_ch)

} else if (params.skip_basecalling && !params.skip_alignment) {

    basecalled_ch = samples_ch
        .filter { it[2] == 'FASTQ' }
        .map { sample, fastq, _ -> tuple(sample, fastq) }

}
    /*
     * STAGE 2
     */
if (!params.skip_alignment) {

    fastq_ch = BAMTOFQ(basecalled_ch)
    aligned_ch = MINIMAP2(fastq_ch)
    sorted_bam_ch = SAMTOOLS_BAM(aligned_ch)
    target_bam_ch = SAMTOOLS_TARGET(sorted_bam_ch)

    if (!params.skip_QC) {
        FLAGSTAT(target_bam_ch)
        CRAMINO(target_bam_ch)
        NANOPLOT(target_bam_ch)
    }

    if (!params.skip_coverage) {
        COVERAGE(target_bam_ch)
    }

} else {

    target_bam_ch = samples_ch
        .filter { it[-1] == 'BAM' }        
        .map { sample, bam, bai, _ -> tuple(sample, bam, bai) }

}
    /*
     * STAGE 3
     */
if (!params.skip_SNVs_PEPPER) {

    snv_results_ch = SNVs_PEPPER(target_bam_ch)
    haplotagged_bam_ch = snv_results_ch[3]   
    INDEX_HP(haplotagged_bam_ch)

}

BIGCLIPPER(target_bam_ch)

if (!params.skip_CALL_SV) {
    CALL_SV(target_bam_ch)
}

}

