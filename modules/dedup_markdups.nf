process DEDUP_MARKDUPS {
  tag "${subject}_${sample_id}"
  container 'quay.io/biocontainers/picard:3.4.0--hdfd78af_0'

  // publish per subject
  publishDir { "${params.outdir_abs}/${subject}/bam" }, mode: 'copy'

  input:
  // match upstream: (subject, sample_id, bam, bai)
  tuple val(subject), val(sample_id), path(bam), path(bai)

  output:
  // keep tuple shape + add subject for downstream joins
  tuple val(subject),
        val(sample_id),
        path("${sample_id}_dedup.bam"),
        path("${sample_id}_dedup.bam.bai"),
        emit: bam

  // metrics for MultiQC etc.
  path "${sample_id}_dedup_metrics.txt", emit: metrics

  script:
  """
  set -euo pipefail

  picard MarkDuplicates \
    I=${bam} \
    O=${sample_id}_dedup.bam \
    M=${sample_id}_dedup_metrics.txt \
    REMOVE_DUPLICATES=true \
    ASSUME_SORTED=true \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT

  # Picard creates ${sample_id}_dedup.bai beside the BAM; normalize to .bam.bai if needed
  if [ -s "${sample_id}_dedup.bai" ] && [ ! -e "${sample_id}_dedup.bam.bai" ]; then
    ln -sf "${sample_id}_dedup.bai" "${sample_id}_dedup.bam.bai"
  fi

  # Sanity checks
  test -s "${sample_id}_dedup.bam"
  test -s "${sample_id}_dedup.bam.bai"
  test -s "${sample_id}_dedup_metrics.txt"
  """
}
