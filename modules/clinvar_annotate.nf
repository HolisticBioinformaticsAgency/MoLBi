process CLINVAR_ANNOTATE {
  tag "${subject}:${id}"
  container '/fs04/scratch2/cg90/liem/projects/nextflow/molbi_pipeline/containers/snpsift_4.3.sif'

  publishDir { "${params.outdir_abs}/${subject}/snpeff" }, mode: 'copy'

  input:
  tuple val(subject), val(id), path(mane_vcf), path(clinvar_vcf)

  output:
  tuple val(subject), val(id), path("${id}_MANE_*.vcf.gz"), emit: snpeffvcf

  script:
  """
  set -euo pipefail

  clin_stem=\$(basename "${clinvar_vcf}")
  clin_stem="\${clin_stem%.gz}"
  clin_stem="\${clin_stem%.vcf}"

  OUT_PREFIX="${id}_MANE_\${clin_stem}"
  OUT_VCF="\${OUT_PREFIX}.vcf"
  OUT_GZ="\${OUT_VCF}.gz"

  # Try fast path first
  set +e
  SnpSift annotate -tabix "${clinvar_vcf}" "${mane_vcf}" > "\${OUT_VCF}"
  rc=\$?
  set -e

  if [[ \$rc -ne 0 ]]; then
    echo "[CLINVAR_ANNOTATE] WARN: -tabix annotate failed. Falling back..." >&2
    gunzip -c "${clinvar_vcf}" > "${id}_clinvar_db.vcf"
    SnpSift annotate -sorted "${id}_clinvar_db.vcf" "${mane_vcf}" > "\${OUT_VCF}"
  fi

  # Compress the final VCF
  gzip -f "\${OUT_VCF}"
  """
}
