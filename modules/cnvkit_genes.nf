process CNVKIT_GENES {
  tag "${subject}_${case_id}"
  container 'quay.io/biocontainers/cnvkit:0.9.12--pyhdfd78af_1'
  env = [ 'MPLCONFIGDIR': '/tmp' ]

  // publish under the case directory
  publishDir { "${pub_base}/cnvkit" }, mode: 'copy'

  input:
  // (pub_base, subject, case_id, cnr, cns)
  tuple val(pub_base), val(subject), val(case_id), path(cnr), path(cns)

  output:
  tuple val(pub_base), val(subject), val(case_id), path("${case_id}_gene_cn.txt"), emit: genecn

  script:
  """
  set -euo pipefail
  cnvkit.py genemetrics \
    ${cnr} \
    -s ${cns} \
    -t ${params.cnvkit_genemetrics_threshold} \
    -m ${params.cnvkit_genemetrics_min_probes} \
    -o ${case_id}_gene_cn.txt
  """
}
