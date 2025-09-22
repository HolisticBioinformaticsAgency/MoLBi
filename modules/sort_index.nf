process SORT_INDEX {
  tag "${subject}_${sample_id}"
  container 'quay.io/biocontainers/samtools:1.18--h50ea8bc_1'

  // publish per subject
  publishDir { "${params.outdir_abs}/${subject}/bam" }, mode: 'copy'

  input:
  tuple val(subject), val(sample_id), path(sam)

  output:
  tuple val(subject), val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam

  script:
  """
  samtools sort -@ ${task.cpus} -o ${sample_id}.bam ${sam}
  samtools index -@ ${task.cpus} ${sample_id}.bam
  """
}
