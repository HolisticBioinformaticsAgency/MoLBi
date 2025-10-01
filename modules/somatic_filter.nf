process SOMATIC_FILTER {
  tag "${subject}:${id}"
  container '/fs04/scratch2/cg90/liem/projects/nextflow/molbi_pipeline/containers/snpsift_4.3.sif'

  publishDir { "${params.outdir_abs}/${subject}/somatic" }, mode: 'copy'

  input:
  tuple val(subject), val(id), path(vcf)

  output:
  tuple val(subject), val(id), path("${id}_somatic_filtered.vcf.gz"), emit: somatic

  script:
  """
  set -euo pipefail

  SnpSift filter "(STATUS = 'StrongLOH') | (STATUS = 'StrongSomatic')" "${vcf}" \
    > "${id}_somatic_filtered.vcf"

  gzip -f "${id}_somatic_filtered.vcf"
  """
}
