nextflow.enable.dsl=2
params.help = false

if (params.help) {
    log.info """
    -----------------------------------------------------------------------
    Basecalling for Long-Read ONT data
    Usage :

nextflow run modules/basecalling.nf -entry basecalling \
    --samplesheet samplesheet.csv \
    --input_pod5 ./test_pod5/ \
    --output_dir . \
    --reference path/to/ref.fa \
    -c nextflow.config \
    --account_name name \
    --use_gpu true \
    --bind_path /path/to/ref/,/path/to/cachedir/,path/to/samples/,etc
 
ONLY WORKS WITH GPU!
    ______________________________________________________________________

    Required:
    
    --samplesheet             CSV file with header:
                               sample,pod5

                              One row per sample, for example:
                               test1,/path/to/pod5/

                              Columns:
                               - sample : sample identifier
                               - pod5   : path to POD5 directory

    --output_dir              Path to output directory
    
    Dorado parameters:

    --dorado_model            Basecalling model with dorado: fast, hac, sup (default "sup")
    --dorado_modified_bases   Space-separated list of modifications following --modified-bases (default "--modified-bases 5mCG_5hmCG,6mA")
    --dorado_params           Other dorado parameters: https://github.com/nanoporetech/dorado/?tab=readme-ov-file (default "--recursive --min-qscore 9 --models-directory /shared/work/PI-tommaso.pippucci/ringtp22/LRS_workflow/dorado_models/")
    --help                    Print this message and exit

       """.stripIndent()
    exit 0
}

if (!params.samplesheet) {
    error "Missing required parameter: --samplesheet"
}

/*
 * Read samplesheet
 */
Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row ->
        def sample = row.sample
        def pod5   = file(row.pod5)

        if (!sample || !pod5)
            error "Invalid row in samplesheet: ${row}"

        tuple(sample, pod5)
    }
    .set { input_pod5 }

params.dorado_model = "sup"
params.dorado_modified_bases = "--modified-bases 5mCG_5hmCG 6mA"
params.dorado_params = "--recursive --min-qscore 9 --models-directory /shared/work/PI-tommaso.pippucci/ringtp22/LRS_workflow/dorado_models/"

process BASECALLING {
    tag "$sample"
    publishDir "${params.output_dir}/results/${sample}", mode: 'copy'
    label (params.use_gpu ? 'gpu' : 'big_job')
   
    input:
    tuple val(sample), path(input_dir)

    output:
    tuple val(sample), path("${sample}.dorado.unaligned.bam")

    script:

def device_flag = params.use_gpu ? '--device cuda:0' : ''

    """
    dorado basecaller ${params.dorado_model} ${input_dir} ${params.dorado_params} ${device_flag} ${params.dorado_modified_bases} > ${sample}.dorado.unaligned.bam
    """
}

process BASECALLING_REPORT {

    tag "${sample}"
    publishDir "${params.output_dir}/results/${sample}", mode: 'copy'
    label 'small_job'

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
	    report = BASECALLING_REPORT(basecalled)
}
