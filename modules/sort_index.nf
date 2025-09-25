process SORT_INDEX {
  tag "${subject}_${sample_id}"
  container 'quay.io/biocontainers/samtools:1.18--h50ea8bc_1'

  // publish per subject
  publishDir { "${params.outdir_abs}/${subject}/bam" }, mode: 'copy'

  input:
  tuple val(subject), val(sample_id), path(sam)

  output:
  tuple val(subject), val(sample_id),
        path("${sample_id}.sorted.bam"),
        path("${sample_id}.sorted.bam.bai"),
        emit: bam

  script:
  """
  set -euo pipefail

  # Keep all reads, no filtering — stream uncompressed BAM to sorter
  samtools view -u -h "${sam}" \
  | samtools sort -@ ${task.cpus} -o "${sample_id}.sorted.bam" -

  # Index the sorted BAM
  samtools index -@ ${task.cpus} "${sample_id}.sorted.bam" "${sample_id}.sorted.bam.bai"

  # Sanity checks
  test -s "${sample_id}.sorted.bam"
  test -s "${sample_id}.sorted.bam.bai"
  """
}
