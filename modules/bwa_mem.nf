process BWA_MEM {
  tag "${subject}_${sample_id}"
  container 'quay.io/biocontainers/bwa:0.7.17--he4a0461_11'

  publishDir { "${params.outdir_abs}/${subject}/bwa_mem" }, mode: 'copy'

  input:
  // tuple = (subject, sample_id, reads, staged FASTA, ABSOLUTE original fasta path string)
  tuple val(subject), val(sample_id), path(reads), path(ref_fa), val(ref_src_abs)

  output:
  tuple val(subject), val(sample_id), path("${sample_id}.sam"), emit: sam

  script:
  """
  set -euo pipefail

  # Symlink BWA sidecar index files next to the staged FASTA so bwa can find them.
  # We use the ABSOLUTE original FASTA path so relative paths can't break.
  for ext in amb ann bwt pac sa; do
    if [ -s "${ref_src_abs}.\$ext" ]; then
      ln -sf "${ref_src_abs}.\$ext" .
    else
      echo "Missing BWA index: ${ref_src_abs}.\$ext" >&2
      exit 2
    fi
  done

  # Read group string (tabs must be escaped for bwa to parse)
  RG="@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPU:lib1\\tPL:ILLUMINA"

  # -M keeps compatibility with Picard/GATK by marking shorter split hits as secondary
  bwa mem -M -t ${task.cpus} -R "\$RG" ${ref_fa} ${reads.join(' ')} > ${sample_id}.sam

  test -s "${sample_id}.sam"
  """
}
