process REF_FAIDX {
  tag { ref_fa.baseName }
  container = 'quay.io/biocontainers/samtools:1.18--h50ea8bc_1'
  publishDir "${params.outdir}/ref_index", mode: 'copy'

  input:
  path ref_fa

  output:
  path "ref.fa",     emit: fasta
  path "ref.fa.fai", emit: fai

  script:
  """
  samtools faidx ${ref_fa}
  """
}

