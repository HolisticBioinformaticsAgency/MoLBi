process SORT_INDEX {
  tag "${subject}_${sample_id}"
  container 'quay.io/biocontainers/samtools:1.18--h50ea8bc_1'

  // publish per subject
  publishDir { "${params.outdir_abs}/${subject}/bam" }, mode: 'copy'

  input:
  tuple val(subject), val(sample_id), path(sam)

  output:
  tuple val(subject), val(sample_id), path("${sample_id}.hq.sorted.bam"), path("${sample_id}.hq.sorted.bam.bai"), emit: bam

  script:
  """
  set -euo pipefail

  # Filter to high-quality, properly paired, mapped primary reads, keep header,
  # stream uncompressed BAM to sorter for speed.
  samtools view -u -h \
      -q 1 \
      -f 2 \
      -F 4 -F 8 -F 256 \
      "${sam}" \
  | samtools sort -@ ${task.cpus} -o "${sample_id}.hq.sorted.bam" -

  # Index the filtered, sorted BAM
  samtools index -@ ${task.cpus} "${sample_id}.hq.sorted.bam" "${sample_id}.hq.sorted.bam.bai"

  # Sanity checks
  test -s "${sample_id}.hq.sorted.bam"
  test -s "${sample_id}.hq.sorted.bam.bai"
  """
}
