process CNVKIT_GENES {
  tag "${subject}_${sample_id}"
  container 'quay.io/biocontainers/cnvkit:0.9.12--pyhdfd78af_1'
  env = [ 'MPLCONFIGDIR': '/tmp' ]

  publishDir { "${params.outdir_abs}/${subject}/cnvkit" }, mode: 'copy'

  input:
  tuple val(subject), val(sample_id), path(cnr), path(cns)

  output:
  tuple val(subject), val(sample_id), path("${sample_id}_gene_cn.txt"), emit: genecn

  script:
  """
  set -euo pipefail
  cnvkit.py genemetrics \
    ${cnr} \
    -s ${cns} \
    -o ${sample_id}_gene_cn.txt
  """
}
