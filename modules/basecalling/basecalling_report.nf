nextflow.enable.dsl=2

process BASECALLING_REPORT {

    tag "${sample}"
    publishDir "${params.output_dir}/results/${sample}", mode: 'copy'
    label 'small_job'

    input:
    tuple val(sample), path(ubam)

    output:
    path("${sample}.basecalling.report.txt")

    script:
    """
set -euo pipefail

UBAM=${ubam}
OUT=${sample}.basecalling.report.txt

    # modello di basecalling
    MODEL=\$(samtools view -H "\$UBAM" | grep -m1 '^@RG' | sed -n 's/.*model=\\([^[:space:]]*\\).*/\\1/p' || true)
    if [ -z "\$MODEL" ]; then
      RG_VAL=\$(samtools view "\$UBAM" | awk 'NR==1{for(i=12;i<=NF;i++) if (\$i ~ /^RG:Z:/) {sub(/^RG:Z:/, "", \$i); print \$i; exit}}' || true)
      if [ -n "\$RG_VAL" ]; then
         MODEL=\$(echo "\$RG_VAL" | sed 's/^[^_]*_//')
      else
         MODEL="NA"
      fi
    fi

    # fallback dal nome file se ancora NA
    if [ "\$MODEL" = "NA" ]; then
      MODEL=\$(basename "\$UBAM" | sed 's/\\.dorado\\.bam//')
    fi

    # modified bases
        MODBASES=\$(samtools view "\$UBAM" \
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
    samtools view -@ ${task.cpus} "\$UBAM" \
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

