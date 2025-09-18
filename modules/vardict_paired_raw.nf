process VARDICT_PAIRED_RAW {
  tag "${case_id}"
  container 'quay.io/biocontainers/vardict-java:1.8.3--hdfd78af_0'
  env     = [ 'JAVA_TOOL_OPTIONS': '-Xms2g -Xmx20g' ]   

  publishDir "${params.outdir}/vardict_paired", mode: 'copy'

  input:
  tuple val(case_id),
        path(t_bam), path(t_bai),
        path(n_bam), path(n_bai),
        path(ref_fa), path(ref_fai),
        path(bed)

  output:
  tuple val(case_id), path("${case_id}.vardict.tsv"), emit: tsv

  script:
  """
  set -euo pipefail
  GCOL=\$(awk 'BEGIN{c=0} !/^#/ && NF{print (NF>=4?4:0); exit}' ${bed})
  echo "[VarDict PAIRED] Using -g \${GCOL}" >&2

  vardict-java \\
    -G ${ref_fa} \\
    -N ${case_id} \\
    -b "${t_bam}|${n_bam}" \\
    -th ${task.cpus} \\
    -f ${params.min_af} -Q 20 -z -c 1 -S 2 -E 3 -g \${GCOL} ${bed} \\
    > ${case_id}.vardict.tsv

  test -s ${case_id}.vardict.tsv || { echo "[VarDict PAIRED] Empty TSV for ${case_id}" >&2; exit 2; }
  """
}
