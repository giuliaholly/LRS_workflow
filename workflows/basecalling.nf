nextflow.enable.dsl=2
params.help = false

if (params.help) {
    log.info """
    -----------------------------------------------------------------------
    Basecalling for Long-Read ONT data
    Usage :
    
 NXF_APPTAINER_CACHEDIR=/shared/work/PI-tommaso.pippucci/ringtp22/my_singularity_container/ NXF_TEMP=/shared/work/PI-tommaso.pippucci/schedulers/tmp/ APPTAINER_TMPDIR=/shared/work/PI-tommaso.pippucci/schedulers/tmp/ nextflow run /shared/work/PI-tommaso.pippucci/ringtp22/basecalling.nf -c /shared/work/PI-tommaso.pippucci/ringtp22/basecalling_nextflow.config --sample sample1 --input /path/to/pod5_directory --output_dir /path/to/output --bind_path
 
ONLY WORKS WITH GPU!
    ______________________________________________________________________

    Required:
    
    --input_pod5              Path to directory of pod5 files 
    --output_dir              Path to output directory
    --sample		      Name of the sample
    
    Dorado parameters:
    --dorado_model            Basecalling model with dorado: fast, hac, sup (default "sup")
    --dorado_modified_bases   Space-separated list of modifications following --modified-bases (default "--modified-bases 5mCG_5hmCG,6mA")
    --dorado_params           Other dorado parameters: https://github.com/nanoporetech/dorado/?tab=readme-ov-file (default "--recursive --min-qscore 9 --models-directory /shared/work/PI-tommaso.pippucci/ringtp22/LRS_workflow/dorado_models/")
    --correct 		      default null
    --help                    Print this message and exit

       """.stripIndent()
    exit 0
}


def samples = params.sample instanceof List ? params.sample : [params.sample]
def inputs = params.input_pod5 instanceof List ? params.input_pod5 : [params.input_pod5]

if (samples.size() != inputs.size()) {
    error "Il numero di samples (${samples.size()}) deve corrispondere al numero di input (${inputs.size()})"
}

def tuples = []
for (int i = 0; i < samples.size(); i++) {
    tuples << [ samples[i], file(inputs[i]) ]
}

Channel.from(tuples).set { input_pod5 }

params.dorado_model = "sup"
params.dorado_modified_bases = "--modified-bases 5mCG_5hmCG 6mA"
params.dorado_params = "--recursive --min-qscore 9 --models-directory ./dorado_models/"
params.correct = null

process BASECALLING {
    tag "$sample"
    publishDir "${params.output_dir}/results/${sample}", mode: 'copy'
    label 'gpu'
   
    input:
    tuple val(sample), path(input_dir)

    output:
    tuple val(sample), path("${sample}.dorado.unaligned.bam")

    script:
    """
	export CUDA_VISIBLE_DEVICES=0
	dorado basecaller ${params.dorado_model} ${input_dir} ${params.dorado_params} --device cuda:0 ${params.dorado_modified_bases} > ${sample}.dorado.unaligned.bam
   
    """
}

process BAMTOFQ {
    tag "$sample"
    publishDir "${params.output_dir}/results/${sample}", mode: 'copy'
    label 'cpu'

    input:
    tuple val(sample), path(sample_input)

    output:
    tuple val(sample), path("${sample}.fastq")

    script:
    """
    samtools fastq -T "*" ${sample_input} > ${sample}.fastq
    
    """
}

process CORRECT {
    tag "$sample"
    publishDir "${params.output_dir}/results/${sample}", mode: 'copy'
    label 'gpu'
  
    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}_corrected.dorado.fasta")

    script:
    """
	export CUDA_VISIBLE_DEVICES=0

	split -l 400000 ${fastq} ${sample}_chunk_

	for f in ${sample}_chunk_*; do	

  		filename=\$(basename "\$f")               
  		output_file="\${filename}_corrected.fasta"

		dorado correct --device cuda:0 "\$f" > ./\"\$output_file\"

	done

	cat \$(ls ${sample}_chunk_*_corrected.fasta | sort) > ${sample}_corrected.dorado.fasta

    """
}

process SAMTOOLS_FAIDX {
    tag "$sample"
    publishDir "${params.output_dir}/results/${sample}", mode: 'copy'
    label 'cpu'

    input:
    tuple val(sample), path("${sample}_corrected.dorado.fasta")

    output:
    tuple val(sample), path("${sample}_corrected.dorado.fasta.fai")
    
    script:
    """

      samtools faidx ${sample}_corrected.dorado.fasta     
    """
}

process BASECALLING_REPORT {

    tag "${sample}"
    publishDir "${params.output_dir}/results/${sample}", mode: 'copy'
    label 'cpu'

    input:
    tuple val(sample), path("${sample}.dorado.unaligned.bam")

    output:
    path("${sample}.basecalling.report.txt")

    script:
    """
set -euo pipefail

BAM=${sample}.dorado.unaligned.bam
OUT=${sample}.basecalling.report.txt

    # modello di basecalling
    MODEL=\$(samtools view -H "\$BAM" | grep -m1 '^@RG' | sed -n 's/.*model=\\([^[:space:]]*\\).*/\\1/p' || true)
    if [ -z "\$MODEL" ]; then
      RG_VAL=\$(samtools view "\$BAM" | awk 'NR==1{for(i=12;i<=NF;i++) if (\$i ~ /^RG:Z:/) {sub(/^RG:Z:/, "", \$i); print \$i; exit}}' || true)
      if [ -n "\$RG_VAL" ]; then
         MODEL=\$(echo "\$RG_VAL" | sed 's/^[^_]*_//')
      else
         MODEL="NA"
      fi
    fi

    # fallback dal nome file se ancora NA
    if [ "\$MODEL" = "NA" ]; then
      MODEL=\$(basename "\$BAM" | sed 's/\\.dorado\\.bam//')
    fi

    # modified bases
        MODBASES=\$(samtools view "\$BAM" \
        | awk '{
            for (i=12; i<=NF; i++) {
                if (\$i ~ /^MM:/ || \$i ~ /^ML:/) {
                    print "ON"; exit
                }
            }
        }
        END { if (!NR) print "OFF" }' \
        || true)

    if [ "\$MODBASES" = "ON" ]; then
      MODBASES="ON"
    else
      MODBASES="OFF"
    fi

    # statistiche
    samtools view -@ ${task.cpus} "\$BAM" \
      | awk -v model="\$MODEL" -v mod="\$MODBASES" '
        BEGIN {
          total_reads=0; total_bases=0; sum_len=0;
          min_len=1e12; max_len=0;
          sum_q=0; q_count=0; min_q=1e12; max_q=-1e12;
        }
        {
          len = length(\$10);
          total_reads++;
          total_bases += len;
          sum_len += len;
          if (len < min_len) min_len = len;
          if (len > max_len) max_len = len;

          q="NA";
          if (match(\$0, /qs:f:-?[0-9]+(\\.[0-9]+)?/)) {
            q = substr(\$0, RSTART+5, RLENGTH-5);
          } else if (match(\$0, /rq:f:-?[0-9]+(\\.[0-9]+)?/)) {
            q = substr(\$0, RSTART+5, RLENGTH-5);
          }

          if (q != "NA") {
            q_val = q+0;
            sum_q += q_val;
            q_count++;
            if (q_val < min_q) min_q=q_val;
            if (q_val > max_q) max_q=q_val;
          }
        }
        END {
          mean_len = (total_reads>0)?sum_len/total_reads:0;
          mean_q = (q_count>0)?sum_q/q_count:"NA";
          if (min_len==1e12) min_len="NA";
          if (min_q==1e12) min_q="NA";
          if (max_q==-1e12) max_q="NA";

          print "Basecall model: " model;
          print "Modified bases: " mod;
          print "Total reads: " total_reads;
          print "Total bases: " total_bases;
          print "Mean read length: " mean_len;
          print "Read length min: " min_len " max: " max_len;
          print "Mean Q-score: " mean_q;
          print "Q-score min: " min_q " max: " max_q;
        }
      ' > "\$OUT"
    """
}

workflow basecalling {
            basecalled = BASECALLING(input_pod5)
            fastq_ready = BAMTOFQ(basecalled)
		if (params.correct) {
	    		fasta = CORRECT(fastq_ready)
			fai = SAMTOOLS_FAIDX(fasta)
   		 } else {
        		log.info "Skipping dorado correct"
    		}

	    report = BASECALLING_REPORT(basecalled)
}
