process PYTMB {
  tag "${subject}:${case_id}"

  container 'quay.io/biocontainers/tmb:1.5.0--pyhdfd78af_1'
  publishDir { "${pub_base}/tmb" }, mode: 'copy'

  input:
    tuple val(pub_base), val(subject), val(case_id), val(tumor_label), path(vcf)
    path db_config
    path var_config

  output:
    tuple val(pub_base), val(subject), val(case_id), path("${case_id}_pytmb.tsv"), emit: tmb

  script:
  """
  set -euo pipefail

  cp ${db_config} db.yml
  cp ${var_config} var.yml

  pyTMB.py \
    -i ${vcf} \
    --effGenomeSize ${params.pytmb_eff_size} \
    --sample "${tumor_label}" \
    --dbConfig db.yml \
    --varConfig var.yml \
    --vaf ${params.pytmb_vaf} \
    --minDepth ${params.pytmb_min_depth} \
    --minAltDepth ${params.pytmb_min_alt} \
    ${params.pytmb_filter_syn ? "--filterSyn" : ""} \
    ${params.pytmb_filter_nc  ? "--filterNonCoding" : ""} \
    ${params.pytmb_filter_lq  ? "--filterLowQual" : ""} \
    > ${case_id}_pytmb.tsv
  """
}