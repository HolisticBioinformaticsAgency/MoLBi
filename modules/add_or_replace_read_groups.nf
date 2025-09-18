process ADD_OR_REPLACE_READ_GROUPS {
  tag "${sample_id}"
  container = 'quay.io/biocontainers/picard:3.4.0--hdfd78af_0'
  publishDir "${params.outdir}/bam", mode: 'copy'

  input:
  tuple val(sample_id), path(bam), path(bai)

  output:
  tuple val(sample_id),
        path("${sample_id}.rg.bam"),
        path("${sample_id}.rg.bam.bai"),
        emit: bam

  script:
  """
  set -euo pipefail

  picard AddOrReplaceReadGroups \
    I=${bam} \
    O=${sample_id}.rg.bam \
    RGID=${sample_id} \
    RGLB=${sample_id} \
    RGPL=ILLUMINA \
    RGPU=${sample_id}.PU \
    RGSM=${sample_id} \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT

  # Picard writes .bai next to bam; normalize to .bam.bai if needed
  if [ -s "${sample_id}.rg.bai" ] && [ ! -e "${sample_id}.rg.bam.bai" ]; then
    ln -sf "${sample_id}.rg.bai" "${sample_id}.rg.bam.bai"
  fi

  test -s "${sample_id}.rg.bam"
  test -s "${sample_id}.rg.bam.bai"
  """
}

