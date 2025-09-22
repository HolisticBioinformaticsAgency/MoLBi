process MSIPRO_MSI {
  tag "${subject} (T:${t_sample}|N:${n_sample})"
  container 'quay.io/biocontainers/msisensor-pro:1.3.0--hd979922_1'

  // per-subject outputs
  publishDir { "${params.outdir_abs}/${subject}/msisensor" }, mode: 'copy'

  input:
  // (subject, t_sample, t_bam, t_bai, n_sample, n_bam, n_bai, sites_tsv, fasta)
  tuple val(subject),
        val(t_sample), path(tumor_bam),  path(tumor_bai),
        val(n_sample), path(normal_bam), path(normal_bai),
        path(sites_tsv),
        path(fasta)

  output:
  // one MSI summary per subject/case
  tuple val(subject), path("${subject}.msi"), emit: msi

  script:
  """
  set -euo pipefail

  # Ensure BAM indexes are discoverable next to BAMs
  [ -e "${tumor_bam}.bai"  ] || ln -sf "${tumor_bai}"  "${tumor_bam}.bai"
  [ -e "${normal_bam}.bai" ] || ln -sf "${normal_bai}" "${normal_bam}.bai"

  # Run msisensor-pro (may output <prefix>, <prefix>.msi, or <prefix>.msi.txt)
  msisensor-pro msi \
    -d ${sites_tsv} \
    -n ${normal_bam} \
    -t ${tumor_bam} \
    -g ${fasta} \
    -o ${subject} \
    -p ${task.cpus} \
    > ${subject}.log 2>&1 || { echo "msisensor-pro failed"; exit 1; }

  # Normalise summary filename to <subject>.msi
  if   [ -s "${subject}.msi" ]; then
    :
  elif [ -s "${subject}.msi.txt" ]; then
    mv "${subject}.msi.txt" "${subject}.msi"
  elif [ -s "${subject}" ] && grep -q '^Total_Number_of_Sites' "${subject}"; then
    mv "${subject}" "${subject}.msi"
  elif [ -s "${subject}.txt" ] && grep -q '^Total_Number_of_Sites' "${subject}.txt"; then
    mv "${subject}.txt" "${subject}.msi"
  else
    echo "ERROR: Could not find msisensor-pro summary; got:" >&2
    ls -l >&2
    echo "==== ${subject}.log (tail) ====" >&2
    tail -n 200 "${subject}.log" >&2 || true
    exit 1
  fi

  test -s "${subject}.msi"
  """
}
