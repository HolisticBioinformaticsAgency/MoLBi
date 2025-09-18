process CLINVAR_ANNOTATE {
  tag "${sample_id}"
  // Your dedicated SnpSift image
  container = '/fs04/scratch2/cg90/liem/projects/nextflow/molbi_pipeline/containers/snpsift_4.3.sif'
  publishDir "${params.outdir}/snpeff", mode: 'copy'

  input:
  // (id, snpeff_core_vcf, clinvar_vcf.gz)
  tuple val(sample_id), path(snpeff_core_vcf), path(clinvar_vcf)

  output:
  tuple val(sample_id), path("${sample_id}.snpeff.vcf"), emit: snpeffvcf

  script:
  """
  set -euo pipefail

  # --- sanity: SnpSift present
  if ! command -v SnpSift >/dev/null 2>&1 ; then
    echo "[CLINVAR_ANNOTATE] ERROR: SnpSift not found in PATH inside container." >&2
    exit 127
  fi

  # --- sanity: ClinVar file looks like a gzip
  if ! gzip -t "${clinvar_vcf}" 2>/dev/null ; then
    echo "[CLINVAR_ANNOTATE] ERROR: ClinVar DB '${clinvar_vcf}' is not a valid .vcf.gz (gzip test failed)." >&2
    exit 2
  fi

  # Try fast path first: use tabix mode on the gz+tabix pair
  set +e
  SnpSift annotate -tabix "${clinvar_vcf}" "${snpeff_core_vcf}" > "${sample_id}.snpeff.vcf"
  rc=\$?
  set -e

  if [[ \$rc -ne 0 ]]; then
    echo "[CLINVAR_ANNOTATE] WARN: -tabix annotate failed (rc=\$rc). Falling back to uncompressed DB..." >&2

    # Fallback: uncompress ClinVar to plain VCF and annotate from that
    # (uses stream-safe gunzip; avoid temp file if very large disk is a concern)
    gunzip -c "${clinvar_vcf}" > "${sample_id}.clinvar.db.vcf"

    # Quick header sanity (should start with '##fileformat=VCF')
    if ! head -n1 "${sample_id}.clinvar.db.vcf" | grep -q '##' ; then
      echo "[CLINVAR_ANNOTATE] ERROR: Uncompressed ClinVar VCF header looks invalid." >&2
      exit 3
    fi

    # If you know ClinVar is sorted (it usually is), hint SnpSift for speed:
    SnpSift annotate -sorted "${sample_id}.clinvar.db.vcf" "${snpeff_core_vcf}" > "${sample_id}.snpeff.vcf"
  fi
  """
}
