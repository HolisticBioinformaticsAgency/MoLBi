process DEDUP_MARKDUPS {
  tag "${sample_id}"
  container = 'quay.io/biocontainers/picard:3.4.0--hdfd78af_0'
  publishDir "${params.outdir}/bam", mode: 'copy'

  input:
  // Same tuple shape as your SORT_INDEX output
  tuple val(sample_id), path(bam), path(bai)

  output:
  // Keep the same tuple shape for downstream tools
  tuple val(sample_id),
        path("${sample_id}.dedup.bam"),
        path("${sample_id}.dedup.bam.bai"),
        emit: bam

  // Also emit metrics for QC dashboards
  path "${sample_id}.dedup.metrics.txt", emit: metrics

  script:
  """
  set -euo pipefail

  picard MarkDuplicates \
    I=${bam} \
    O=${sample_id}.dedup.bam \
    M=${sample_id}.dedup.metrics.txt \
    REMOVE_DUPLICATES=true \
    ASSUME_SORTED=true \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT

  # Picard creates ${sample_id}.dedup.bai beside the BAM; normalize to .bam.bai
  if [ -s "${sample_id}.dedup.bai" ] && [ ! -e "${sample_id}.dedup.bam.bai" ]; then
    ln -sf "${sample_id}.dedup.bai" "${sample_id}.dedup.bam.bai"
  fi

  # Sanity
  test -s "${sample_id}.dedup.bam"
  test -s "${sample_id}.dedup.bam.bai"
  test -s "${sample_id}.dedup.metrics.txt"
  """
}

