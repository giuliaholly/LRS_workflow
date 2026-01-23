nextflow.enable.dsl = 2

if (params.help) {
    log.info "Use --samplesheet samplesheet.tsv"
    exit 0
}

if (!params.samplesheet) {
    error "You must provide --samplesheet"
}

params.bed = false

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
    
    bam_for_qc_ch = params.bed 
        ? SAMTOOLS_TARGET(sorted_bam_ch) 
        : sorted_bam_ch
 
    if (!params.skip_QC) {
        FLAGSTAT(bam_for_qc_ch)
        CRAMINO(bam_for_qc_ch)
        NANOPLOT(bam_for_qc_ch)
    }

    if (!params.skip_coverage) {
        COVERAGE(bam_for_qc_ch)
    }

} else {

    bam_for_qc_ch = samples_ch
        .filter { it[-1] == 'BAM' }        
        .map { sample, bam, bai, _ -> tuple(sample, bam, bai) }

}
    /*
     * STAGE 3
     */
    if (!params.skip_SNVs_PEPPER) {
        snv_results_ch = SNVs_PEPPER(bam_for_qc_ch)
        haplotagged_bam_ch = snv_results_ch[3]   
        INDEX_HP(haplotagged_bam_ch)
    }

    BIGCLIPPER(bam_for_qc_ch)

    if (!params.skip_CALL_SV) {
        CALL_SV(bam_for_qc_ch)
    }

}

