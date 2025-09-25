process CNVKIT_GENES {
  tag "${subject}_${sample_id}"
  container 'quay.io/biocontainers/cnvkit:0.9.12--pyhdfd78af_1'
  env = [ 'MPLCONFIGDIR': '/tmp' ]

<<<<<<< HEAD
  // per-subject outputs
=======
>>>>>>> f12105e (Pipeline moved to vh83, dropped filtering flags for sort_index)
  publishDir { "${params.outdir_abs}/${subject}/cnvkit" }, mode: 'copy'

  input:
  tuple val(subject), val(sample_id), path(cnr), path(cns)

  output:
  tuple val(subject), val(sample_id), path("${sample_id}_gene_cn.txt"), emit: genecn

  script:
  """
  set -euo pipefail
<<<<<<< HEAD

=======
>>>>>>> f12105e (Pipeline moved to vh83, dropped filtering flags for sort_index)
  cnvkit.py genemetrics \
    ${cnr} \
    -s ${cns} \
    -o ${sample_id}_gene_cn.txt
  """
}
