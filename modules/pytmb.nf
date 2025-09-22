process PYTMB {
  tag "${subject}:${sample_id}"

  container 'quay.io/biocontainers/tmb:1.5.0--pyhdfd78af_1'
  publishDir "${params.outdir_abs}/${subject}/tmb", mode: 'copy'

  input:
    tuple val(subject), val(sample_id), path(vcf)
    path db_config
    path var_config

  output:
    tuple val(subject), val(sample_id), path("${sample_id}_pytmb.tsv"), emit: tmb

  script:
  """
  set -euo pipefail

  cp ${db_config} db.yml
  cp ${var_config} var.yml

  # Always use tumour column (with "_T") for pyTMB
  pyTMB.py \
    -i ${vcf} \
    --effGenomeSize ${params.pytmb_eff_size} \
    --sample ${sample_id}_T \
    --dbConfig db.yml \
    --varConfig var.yml \
    --vaf ${params.pytmb_vaf} \
    --minDepth ${params.pytmb_min_depth} \
    --minAltDepth ${params.pytmb_min_alt} \
    ${params.pytmb_filter_syn ? "--filterSyn" : ""} \
    ${params.pytmb_filter_nc  ? "--filterNonCoding" : ""} \
    ${params.pytmb_filter_lq  ? "--filterLowQual" : ""} \
    > ${sample_id}_pytmb.tsv
  """
}
