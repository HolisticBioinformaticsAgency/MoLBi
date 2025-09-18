process SNPEFF_ANNOTATE {
  tag "${sample_id}"
  container = 'quay.io/biocontainers/snpeff:5.2--hdfd78af_1'
  publishDir "${params.outdir}/snpeff", mode: 'copy'
  // profile binds /snpeff_data: withName: 'SNPEFF_ANNOTATE' { containerOptions = "--bind /home/.../snpEff:/snpeff_data" }

  input:
  // (id, raw_vcf, MANE_dir)
  tuple val(sample_id), path(vcf), path(mane_db_dir)

  output:
  // Emit an intermediate SnpEff-only VCF (no ClinVar yet)
  tuple val(sample_id), path("${sample_id}.snpeff.core.vcf"), emit: core

  script:
  """
  set -euo pipefail

  BASE_CFG="/snpeff_data/snpEff.config"
  CFG="snpeff.mane.config"

  # start from base cfg if available
  if [[ -s "\$BASE_CFG" ]]; then
    cp "\$BASE_CFG" "\$CFG"
  else
    echo "data.dir = \$PWD" > "\$CFG"
  fi
  # force data.dir to the work dir and register the MANE genome name
  {
    echo ""
    echo "data.dir = \$PWD"
    echo "GRCh38.mane.1.2.refseq.genome : Homo sapiens (GRCh38 MANE v1.2 RefSeq)"
  } >> "\$CFG"

  # ensure SnpEff sees the genome at ./GRCh38.mane.1.2.refseq
  [[ -d "GRCh38.mane.1.2.refseq" ]] || ln -s "${mane_db_dir}" "GRCh38.mane.1.2.refseq"

  # drop symbolic SV alleles (<DEL>, <DUP>, <INV>, <INS>, <CNV>, etc.) that crash SnpEff
  awk 'BEGIN{OFS="\\t"} /^#/ {print; next} \$5 ~ /^</ {next} {print}' "${vcf}" > "${sample_id}.snpeff.input.vcf"

  # run SnpEff on SNVs/indels only using MANE db
  snpEff ann \\
    -c "\$CFG" \\
    -dataDir "\$PWD" \\
    -noStats -no-downstream -no-intergenic -no-upstream \\
    "GRCh38.mane.1.2.refseq" \\
    "${sample_id}.snpeff.input.vcf" > "${sample_id}.snpeff.core.vcf"
  """
}
