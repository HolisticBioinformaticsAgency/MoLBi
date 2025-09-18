process POLYSOLVER {
  tag "${sample_id}"
  container = 'shub://IARCbioinfo/polysolver-singularity:v4'
  publishDir "${params.outdir}/polysolver", mode: 'copy'

  // Ensure inputs are real files in the task dir (not symlinks)
  stageInMode 'copy'

  input:
  tuple val(sample_id), path(bam), path(bai)

  output:
  tuple val(sample_id), path("${sample_id}.hla.txt"), emit: hla

  // Tunables
  def RACE       = params.polysolver_race       ?: 'Caucasian'
  def BUILD      = params.polysolver_build      ?: 'hg38'
  def EMITVCF    = params.polysolver_emit_vcf   ?: 0
  def FQTYPE     = params.polysolver_fastqtype  ?: 'STDFQ'
  def INSERTCALC = params.polysolver_insertcalc ?: 0

  script:
  """
  set -euo pipefail

  WRAPPER="/home/polysolver/scripts/shell_call_hla_type"
  [ -x "\$WRAPPER" ] || { echo "ERROR: \$WRAPPER not found in container"; exit 127; }

  echo "DEBUG: task dir is: \$(pwd)"
  echo "DEBUG: BAM path: ${bam}"
  ls -l "${bam}" || true
  [ -s "${bam}" ] || { echo "ERROR: BAM not found or empty: ${bam}"; exit 1; }

  # Make sure a .bam.bai exists alongside the BAM if the index was staged separately
  if [ ! -e "${bam}.bai" ] && [ -s "${bai}" ]; then
    ln -sf "${bai}" "${bam}.bai" || true
  fi

  # Polysolver positional args:
  #   <bam> <race> <include_vcf> <build> <fastq> <insertCalc> <outDir>
  "\$WRAPPER" \
    "${bam}" \
    "${RACE}" \
    "${EMITVCF}" \
    "${BUILD}" \
    "${FQTYPE}" \
    "${INSERTCALC}" \
    "\$(pwd)" > "${sample_id}.hla.txt"
  """
}
