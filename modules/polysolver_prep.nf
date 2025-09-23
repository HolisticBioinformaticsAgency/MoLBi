process POLYSOLVER_PREP {
  tag "${subject}_${sample_id}"
  container 'quay.io/biocontainers/samtools:1.18--h50ea8bc_1'
  publishDir { "${params.outdir_abs}/${subject}/polysolver" }, mode: 'copy'
  cpus 4

  input:
  tuple val(subject), val(sample_id), path(bam), path(bai), path(fai)

  output:
  tuple val(subject), val(sample_id),
        path("${sample_id}.ps.bam"),
        path("${sample_id}.ps.bam.bai"),
        emit: ps_bam

  script:
  """
  set -euo pipefail
  export LC_ALL=C

  in_bam="${bam}"
  in_bai="${bai}"
  REF_FAI="${fai}"
  tcpus="${task.cpus}"

  # quick guard
  if [ ! -s "\$REF_FAI" ]; then
    echo "ERROR: ref_fai not provided or empty: '\$REF_FAI'" >&2
    ls -l || true
    exit 2
  fi

  # Ensure the index sits next to the BAM (helps idxstats)
  if [ ! -e "\${in_bam}.bai" ] && [ -s "\${in_bai}" ]; then
    ln -sfn "\${in_bai}" "\${in_bam}.bai"
  fi

  # ---- Canonical contigs from FAI (1..22, X, Y, M/MT; +/- 'chr') ----
  awk -F'\\t' 'BEGIN{OFS="\\t"}
    \$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)\$/ { print \$1, \$2 }
  ' "\$REF_FAI" > fai.canon.tsv

  [ -s fai.canon.tsv ] || { echo "ERROR: No canonical contigs matched in FAI: \$REF_FAI" >&2; exit 3; }

  # ---- Intersect with BAM contigs ----
  samtools idxstats "\${in_bam}" | awk -F'\\t' '\$1!="*"{print \$1}' | sort -u > bam.contigs
  cut -f1 fai.canon.tsv | sort -u > fai.contigs
  comm -12 bam.contigs fai.contigs > keep.contigs
  [ -s keep.contigs ] || { echo "ERROR: No overlapping contigs between BAM and FAI" >&2; exit 4; }

  # ---- Minimal header: @HD + @SQ ----
  printf "@HD\\tVN:1.0\\tSO:coordinate\\n" > header.min.sam
  awk -F'\\t' 'NR==FNR{keep[\$1]=1; next} keep[\$1]{printf "@SQ\\tSN:%s\\tLN:%s\\n", \$1, \$2}' \
    keep.contigs fai.canon.tsv >> header.min.sam

  # ---- Subset & reheader ----
  samtools view -@ "\${tcpus}" -b "\${in_bam}" \$(tr '\\n' ' ' < keep.contigs) > body.canon.bam
  samtools reheader -P header.min.sam body.canon.bam > "${sample_id}.ps.bam"
  samtools index -@ "\${tcpus}" "${sample_id}.ps.bam" "${sample_id}.ps.bam.bai"

  test -s "${sample_id}.ps.bam"
  test -s "${sample_id}.ps.bam.bai"
  """
}
