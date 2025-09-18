process VARDICT_TO_VCF {
  tag "${id}"
  container 'quay.io/biocontainers/vardict-java:1.8.3--hdfd78af_0'
  publishDir "${params.outdir}/variants", mode: 'copy'

  input:
  // id = case/sample id; mode = 'single' or 'paired'
  tuple val(id), val(mode), path(tsv)

  output:
  tuple val(id), path("${id}.vcf"), emit: vcf

  script:
  """
  set -euo pipefail
  echo "[VarDict TO_VCF] Mode: ${mode}" >&2

  if [ "${mode}" = "paired" ]; then
    cat ${tsv} \
      | testsomatic.R \
      | var2vcf_paired.pl -N "${id}|${id}_NORMAL" -f ${params.min_af} \
      > ${id}.vcf
  else
    cat ${tsv} \
      | teststrandbias.R \
      | var2vcf_valid.pl -N ${id} -E -f ${params.min_af} \
      > ${id}.vcf
  fi
  """
}
