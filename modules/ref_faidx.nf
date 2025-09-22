process REF_FAIDX {
  tag "ref_faidx"
  container 'quay.io/biocontainers/samtools:1.18--h50ea8bc_1'

  input:
  path ref_fa

  output:
  path "${ref_fa}.fai", emit: fai
  path ref_fa,          emit: fasta

  script:
  """
  set -euo pipefail
  samtools faidx ${ref_fa}
  """
}
