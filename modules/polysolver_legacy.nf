process POLYSOLVER {
  tag "${subject}_${sample_id}"
  container 'shub://IARCbioinfo/polysolver-singularity:v4'
  publishDir { "${params.outdir_abs}/${subject}/polysolver" }, mode: 'copy'

  input:
  tuple val(subject), val(sample_id), path(ps_bam), path(ps_bai)

  output:
  tuple val(subject), val(sample_id), path("${sample_id}_hla.txt"), emit: hla

  script:
  """
  set -euo pipefail
  export TMPDIR="\$(pwd)"
  export JAVA_TOOL_OPTIONS="\${JAVA_TOOL_OPTIONS:-} -Djava.io.tmpdir=\$TMPDIR"
  export LC_ALL=C

  WRAPPER="/home/polysolver/scripts/shell_call_hla_type"
  [ -x "\$WRAPPER" ] || { echo "ERROR: wrapper not found"; exit 127; }

  # Use the staged-in names as provided by Nextflow
  bam="${ps_bam}"
  bai="${ps_bai}"

  echo "DEBUG: BAM=\$bam"
  ls -l "\$bam" || true
  [ -s "\$bam" ] || { echo "ERROR: BAM missing: \$bam"; exit 1; }

  # Ensure an index with the SAME BASENAME sits next to the BAM
  # If the index file isn't already exactly <bam>.bai, create a symlink.
  if [ ! -e "\${bam}.bai" ]; then
    if [ -s "\$bai" ]; then
      ln -s "\$bai" "\${bam}.bai"
    elif [ -s "\${bam}.bai" ]; then
      :  # already there
    else
      echo "ERROR: No usable BAI for \$bam" >&2
      exit 1
    fi
  fi

  [ -s "\${bam}.bai" ] || { echo "ERROR: BAM index missing for \$bam"; exit 1; }

  # Run POLYSOLVER with your chosen parameters (unchanged)
  "\$WRAPPER" \
    "\$bam" \
    "${params.polysolver_race ?: 'Caucasian'}" \
    "${params.polysolver_emit_vcf ?: 0}" \
    "${params.polysolver_build ?: 'hg38'}" \
    "${params.polysolver_fastqtype ?: 'STDFQ'}" \
    "${params.polysolver_insertcalc ?: 0}" \
    "\$(pwd)" > "${sample_id}_hla.txt"

  test -s "${sample_id}_hla.txt"
  """
}
