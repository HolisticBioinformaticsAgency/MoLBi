process SOMATIC_FILTER {
  
  tag "${subject}:${id}"
  container '/fs04/scratch2/cg90/liem/projects/nextflow/molbi_pipeline/containers/snpsift_4.3.sif'

  // let the bgzip wrapper publish later, so no publishDir here (or keep yours if you want)
  // publishDir { "${pub_base}/somatic" }, mode: 'copy'

  input:
  // (pub_base, subject, id, mode, vcf)  — vcf is now very often .vcf.gz
  tuple val(pub_base), val(subject), val(id), val(mode), path(vcf)

  output:
  // always a plain VCF here, bgzip step will handle compression
  tuple val(pub_base), val(subject), val(id), path("${id}_somatic_filtered.vcf"), emit: somatic

  script:
  """
  set -euo pipefail

  # always normalize to plain VCF first
  zcat -f "${vcf}" > "${id}_somatic_input.vcf"

  if [ "${mode}" = "paired" ]; then
    echo "[SomaticFilter] Paired mode: filtering STATUS" >&2
    # feed the *plain* VCF to SnpSift
    SnpSift filter "(STATUS = 'StrongLOH') | (STATUS = 'StrongSomatic')" "${id}_somatic_input.vcf" \
      > "${id}_somatic_filtered.vcf"
  else
    echo "[SomaticFilter] Single mode: skipping STATUS filter" >&2
    cp "${id}_somatic_input.vcf" "${id}_somatic_filtered.vcf"
  fi

  # sanity
  if ! grep -q '^#CHROM' "${id}_somatic_filtered.vcf"; then
    echo "[SomaticFilter] ERROR: output VCF missing #CHROM header" >&2
    head -50 "${id}_somatic_filtered.vcf" >&2
    exit 1
  fi
  """
}
