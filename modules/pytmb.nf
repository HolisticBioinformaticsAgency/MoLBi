process PYTMB {
  tag "${sample_id}"
  container 'quay.io/biocontainers/tmb:1.5.0--pyhdfd78af_1'
  publishDir "${params.outdir}/tmb", mode: 'copy'

  input:
  tuple val(sample_id), path(vcf)
  path db_config
  path var_config

  output:
  tuple val(sample_id), path("${sample_id}.pytmb.tsv"), emit: tmb

  script:
  """
  set -euo pipefail

  # Copy configs into the task dir
  cp ${db_config} db.yml
  cp ${var_config} var.yml

  pyTMB.py \
    -i ${vcf} \
    --effGenomeSize ${params.pytmb_eff_size} \
    --sample ${sample_id} \
    --dbConfig db.yml \
    --varConfig var.yml \
    --vaf ${params.pytmb_vaf} \
    --minDepth ${params.pytmb_min_depth} \
    --minAltDepth ${params.pytmb_min_alt} \
    ${params.pytmb_filter_syn ? "--filterSyn" : ""} \
    ${params.pytmb_filter_nc ? "--filterNonCoding" : ""} \
    ${params.pytmb_filter_lq ? "--filterLowQual" : ""} \
    > ${sample_id}.pytmb.tsv
  """
}
