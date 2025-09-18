process BWA_INDEX {
  tag { ref_fa.baseName }
  container = 'quay.io/biocontainers/bwa:0.7.17--he4a0461_11'
  publishDir "${params.outdir}/ref_index", mode: 'copy'

  input:
  path ref_fa

  output:
  path "ref_index", emit: idxdir
  path "ref.fa",    emit: fasta

  script:
  """
  ln -s ${ref_fa} ref.fa
  mkdir -p ref_index
  bwa index -p ref_index/genome ref.fa
  """
}

