process ADD_OR_REPLACE_READ_GROUPS {
  tag "${subject}_${sample_id}"
  container 'quay.io/biocontainers/picard:3.4.0--hdfd78af_0'

  // publish per subject
  publishDir { "${params.outdir_abs}/${subject}/bam" }, mode: 'copy'

  input:
  tuple val(subject), val(sample_id), path(bam), path(bai)

  output:
  tuple val(subject),
        val(sample_id),
        path("${sample_id}_rg.bam"),
        path("${sample_id}_rg.bam.bai"),
        emit: bam

  script:
  """
  set -euo pipefail

  picard AddOrReplaceReadGroups \
    I=${bam} \
    O=${sample_id}_rg.bam \
    RGID=${sample_id} \
    RGLB=${sample_id} \
    RGPL=ILLUMINA \
    RGPU=${sample_id}.PU \
    RGSM=${sample_id} \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT

  # Picard writes .bai next to BAM; normalize to .bam.bai if needed
  if [ -s "${sample_id}_rg.bai" ] && [ ! -e "${sample_id}_rg.bam.bai" ]; then
    ln -sf "${sample_id}_rg.bai" "${sample_id}_rg.bam.bai"
  fi

  test -s "${sample_id}_rg.bam"
  test -s "${sample_id}_rg.bam.bai"
  """
}
