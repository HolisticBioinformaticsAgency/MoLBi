process SNPEFF_ANNOTATE {
  tag "${subject}:${id}"
  container 'quay.io/biocontainers/snpeff:5.2--hdfd78af_1'

  publishDir { "${params.outdir_abs}/${subject}/snpeff" }, mode: 'copy'

  input:
  tuple val(subject), val(id), path(vcf), path(mane_db_dir)

  output:
  tuple val(subject), val(id), path("${id}_snpeff_core.vcf"), emit: core

  script:
  """
  set -euo pipefail

  BASE_CFG="/snpeff_data/snpEff.config"
  CFG="snpeff.mane.config"

  if [[ -s "\$BASE_CFG" ]]; then
    cp "\$BASE_CFG" "\$CFG"
  else
    echo "data.dir = \$PWD" > "\$CFG"
  fi

  {
    echo ""
    echo "data.dir = \$PWD"
    echo "GRCh38.mane.1.2.refseq.genome : Homo sapiens (GRCh38 MANE v1.2 RefSeq)"
  } >> "\$CFG"

  [[ -d "GRCh38.mane.1.2.refseq" ]] || ln -s "${mane_db_dir}" "GRCh38.mane.1.2.refseq"

  awk 'BEGIN{OFS="\\t"} /^#/ {print; next} \$5 ~ /^</ {next} {print}' "${vcf}" > "${id}_snpeff_input.vcf"

  snpEff ann \\
    -c "\$CFG" \\
    -dataDir "\$PWD" \\
    -noStats -no-downstream -no-intergenic -no-upstream \\
    "GRCh38.mane.1.2.refseq" \\
    "${id}_snpeff_input.vcf" > "${id}_snpeff_core.vcf"
  """
}
