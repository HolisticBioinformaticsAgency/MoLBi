process ALIGN_SAM {
  tag "${sample_id}"
  container = 'quay.io/biocontainers/bwa:0.7.17--he4a0461_11'
  publishDir "${params.outdir}/align_sam", mode: 'copy'

  input:
  tuple val(sample_id), path(reads), path(idxdir), path(ref_fa)

  output:
  tuple val(sample_id), path("${sample_id}.sam"), emit: sam

  script:
  """
  bwa mem -t ${task.cpus} ${idxdir}/genome ${reads.join(" ")} > ${sample_id}.sam
  """
}

