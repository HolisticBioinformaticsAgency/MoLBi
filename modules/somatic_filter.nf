process SOMATIC_FILTER {
  tag "${subject}:${id}"
  container '/fs04/scratch2/cg90/liem/projects/nextflow/molbi_pipeline/containers/snpsift_4.3.sif'

  // per-subject outputs
  publishDir { "${params.outdir_abs}/${subject}/somatic" }, mode: 'copy'

  input:
  // (subject, id, vcf)
  tuple val(subject), val(id), path(vcf)

  output:
  tuple val(subject), val(id), path("${id}_somatic.vcf"), emit: somatic

  script:
  """
  set -euo pipefail

  SnpSift filter "(STATUS = 'StrongLOH') | (STATUS = 'StrongSomatic')" "${vcf}" \
    > "${id}_somatic.vcf"
  """
}
