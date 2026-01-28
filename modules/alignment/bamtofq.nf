nextflow.enable.dsl=2

process BAMTOFQ {
    tag "$sample"
    label 'big_job'

    input:
    tuple val(sample), path(input_file)

    output:
    tuple val(sample), path("${sample}.fastq")

    script:
    def ext = input_file.getName().tokenize('.')[-1]

    if (ext == 'bam') {
        """
        samtools fastq -T "*" ${input_file} > ${sample}.fastq
        """
    } else if (ext == 'fastq') {
        """
        cp ${input_file} ${sample}.fastq
        """
    } else {
        error "Unsupported file extension: $ext"
    }
}

