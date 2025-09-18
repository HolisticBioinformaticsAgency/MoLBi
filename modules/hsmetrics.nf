process HSMETRICS {
  tag "${sample_id}"
  container = 'quay.io/biocontainers/picard:3.4.0--hdfd78af_0'
  publishDir "${params.outdir}/hsmetrics", mode: 'copy'

  input:
  // IMPORTANT: take the BAM from DEDUP_MARKDUPS, not from SORT_INDEX
  tuple val(sample_id), path(bam), path(bai)     // from DEDUP_MARKDUPS.emit: bam
  path bed                                       // e.g. S33266436_Padded.hg38matched.bed
  path ref_fa                                    // FASTA path
  path ref_fai                                   // FASTA .fai (unused by Picard but fine to pass)

  output:
  tuple val(sample_id), path("${sample_id}.hs_metrics.txt"), emit: hs

  script:
  // Use variable names consistently; avoid hard-coded file names
  def dict = "${ref_fa.baseName}.dict"

  """
  set -euo pipefail

  # Build dictionary next to the staged reference if needed
  if [ ! -s "${dict}" ]; then
    picard CreateSequenceDictionary \
      -R ${ref_fa} \
      -O ${dict}
  fi

  # BED -> interval_list matching the reference dictionary
  picard BedToIntervalList \
    -I ${bed} \
    -O targets.interval_list \
    -SD ${dict}

  # Collect HsMetrics on the staged BAM
  picard CollectHsMetrics \
    -I ${bam} \
    -O ${sample_id}.hs_metrics.txt \
    -R ${ref_fa} \
    -BAIT_INTERVALS targets.interval_list \
    -TARGET_INTERVALS targets.interval_list \
    -VALIDATION_STRINGENCY SILENT

  # sanity
  test -s ${sample_id}.hs_metrics.txt
  """
}
