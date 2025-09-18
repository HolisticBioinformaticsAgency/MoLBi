process CNVKIT_GENES {
  tag "${sample_id}"
  container 'quay.io/biocontainers/cnvkit:0.9.12--pyhdfd78af_1'
  env = [ 'MPLCONFIGDIR': '/tmp' ]

  publishDir "${params.outdir}/cnvkit", mode: 'copy'

  input:
  tuple val(sample_id), path(cnr), path(cns)

  output:
  tuple val(sample_id), path("${sample_id}.gene_cn.txt"), emit: genecn

  script:
  """
  cnvkit.py genemetrics \
    ${cnr} \
    -s ${cns} \
    -o ${sample_id}.gene_cn.txt
  """
}
