process POLYSOLVER {
  tag "${subject}_${sample_id}"
  container 'shub://IARCbioinfo/polysolver-singularity:v4'
  publishDir { "${params.outdir_abs}/${subject}/polysolver" }, mode: 'copy'

  input:
  // matches ch_bam_for_hla_in: (subject, sample, bam, bai, ref_fai)
  tuple val(subject), val(sample_id), path(bam), path(bai) 

  output:
  tuple val(subject), val(sample_id), path("${sample_id}_hla.txt"), emit: hla

  script:
  """
  set -euo pipefail
  export LC_ALL=C

  # stay in task dir for all temp
  export TMPDIR="\$(pwd)"
  mkdir -p "\$TMPDIR"
  export JAVA_TOOL_OPTIONS="\${JAVA_TOOL_OPTIONS:-} -Djava.io.tmpdir=\$TMPDIR"

  WRAPPER="/home/polysolver/scripts/shell_call_hla_type"
  [ -x "\$WRAPPER" ] || { echo "ERROR: wrapper not found"; exit 127; }

  echo "DEBUG: task dir: \$(pwd)"
  echo "DEBUG: bam=${bam}"
  ls -l "${bam}" || true
  [ -s "${bam}" ] || { echo "ERROR: BAM missing: ${bam}"; exit 1; }

  # Make sure an index with the SAME BASENAME sits next to the BAM
  if [ ! -e "${bam}.bai" ]; then
    if [ -s "${bai}" ]; then
      ln -s "${bai}" "${bam}.bai"
    elif [ -s "${bam}.bai" ]; then
      :  # already there
    else
      echo "ERROR: No usable BAI for ${bam}" >&2
      exit 1
    fi
  fi
  [ -s "${bam}.bai" ] || { echo "ERROR: BAM index missing for ${bam}"; exit 1; }

  # Run polysolver
  "\$WRAPPER" \
    "${bam}" \
    "${params.polysolver_race       ?: 'Caucasian'}" \
    "${params.polysolver_emit_vcf   ?: 0}" \
    "${params.polysolver_build      ?: 'hg38'}" \
    "${params.polysolver_fastqtype  ?: 'STDFQ'}" \
    "${params.polysolver_insertcalc ?: 0}" \
    "\$(pwd)" > "${sample_id}_hla.txt"

  test -s "${sample_id}_hla.txt"
  """
}
