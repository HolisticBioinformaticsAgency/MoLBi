process SOMATIC_FILTER {
  tag "${sample_id}"
  container = '/fs04/scratch2/cg90/liem/projects/nextflow/molbi_pipeline/containers/snpsift_4.3.sif'
  publishDir "${params.outdir}/somatic", mode: 'copy'

  input:
  tuple val(sample_id), path(vcf)

  output:
  tuple val(sample_id), path("${sample_id}.somatic.vcf"), emit: somatic

  script:
  """
  set -euo pipefail
  SnpSift filter "(STATUS = 'StrongLOH') | (STATUS = 'StrongSomatic')" "${vcf}" > "${sample_id}.somatic.vcf"
  """
}