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
  
  IN_VCF="${vcf}"
  PLAIN_VCF="${case_id}_somatic_filtered.vcf"
  
  # Decompress if needed
  if [[ "\${IN_VCF}" == *.gz ]]; then
    gunzip -c "\${IN_VCF}" > "\${PLAIN_VCF}"
  else
    cp "\${IN_VCF}" "\${PLAIN_VCF}"
  fi
  
  # Validate VCF has content
  if [ ! -s "\${PLAIN_VCF}" ]; then
    echo "ERROR: VCF file is empty" >&2
    exit 1
  fi
  
  # Check if VCF has any header lines (##) - more lenient check
  if ! grep -q "^##" "\${PLAIN_VCF}"; then
    echo "ERROR: VCF file has no header lines" >&2
    head -20 "\${PLAIN_VCF}" >&2
    exit 1
  fi
  
  # Check if VCF has the column header line
  if ! grep -q "^#CHROM" "\${PLAIN_VCF}"; then
    echo "ERROR: VCF file missing #CHROM header line" >&2
    head -20 "\${PLAIN_VCF}" >&2
    exit 1
  fi
  
  # Count variants (excluding headers)
  VARIANT_COUNT=\$(grep -v "^#" "\${PLAIN_VCF}" | wc -l)
  echo "INFO: Found \${VARIANT_COUNT} variants in VCF" >&2
  
  if [ "\${VARIANT_COUNT}" -eq 0 ]; then
    echo "WARNING: No variants found in VCF, creating empty TMB output" >&2
    echo -e "sample\\ttmb\\ttotal_variants\\tcoding_variants" > ${case_id}_pytmb.tsv
    echo -e "${tumor_label}\\t0\\t0\\t0" >> ${case_id}_pytmb.tsv
  else
    pyTMB.py \\
      -i "\${PLAIN_VCF}" \\
      --effGenomeSize ${params.pytmb_eff_size} \\
      --sample "${tumor_label}" \\
      --dbConfig db.yml \\
      --varConfig var.yml \\
      --vaf ${params.pytmb_vaf} \\
      --minDepth ${params.pytmb_min_depth} \\
      --minAltDepth ${params.pytmb_min_alt} \\
      ${params.pytmb_filter_syn ? "--filterSyn" : ""} \\
      ${params.pytmb_filter_nc  ? "--filterNonCoding" : ""} \\
      ${params.pytmb_filter_lq  ? "--filterLowQual" : ""} \\
      > ${case_id}_pytmb.tsv
  fi
  """
}