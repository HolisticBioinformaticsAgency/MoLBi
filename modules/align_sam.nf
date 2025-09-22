process ALIGN_SAM {
  tag "${subject}_${sample_id}"
  container 'quay.io/biocontainers/bwa:0.7.17--he4a0461_11'

  // publish per subject
  publishDir { "${params.outdir_abs}/${subject}/align_sam" }, mode: 'copy'

  input:
  tuple val(subject), val(sample_id), path(reads), path(idxdir), path(ref_fa)

  output:
  tuple val(subject), val(sample_id), path("${sample_id}.sam"), emit: sam

  script:
  """
  bwa mem -t ${task.cpus} ${idxdir}/genome ${reads.join(" ")} > ${sample_id}.sam
  """
}
