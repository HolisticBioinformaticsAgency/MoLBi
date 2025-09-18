process MSIPRO_MSI {
  tag "${sample_id}"
  container = 'quay.io/biocontainers/msisensor-pro:1.3.0--hd979922_1'
  publishDir "${params.outdir}/msisensor", mode: 'copy'

  input:
  tuple val(sample_id),
        path(tumor_bam),  path(tumor_bai),
        path(normal_bam), path(normal_bai),
        path(sites_tsv),
        path(fasta)

  output:
  // ensure Nextflow always finds "<sample_id>.msi"
  tuple val(sample_id), path("${sample_id}.msi"), emit: msi

  script:
  """
  set -euo pipefail

  # Make sure index files are discoverable next to BAMs (name.bam.bai)
  [ -e "${tumor_bam}.bai"  ] || ln -sf "${tumor_bai}"  "${tumor_bam}.bai"
  [ -e "${normal_bam}.bai" ] || ln -sf "${normal_bai}" "${normal_bam}.bai"

  # Run msisensor-pro — note: it may write either "<prefix>" or "<prefix>.msi"
  msisensor-pro msi \
    -d ${sites_tsv} \
    -n ${normal_bam} \
    -t ${tumor_bam} \
    -g ${fasta} \
    -o ${sample_id} \
    -p ${task.cpus} \
    > ${sample_id}.log 2>&1 || { echo "msisensor-pro failed"; exit 1; }

  # Normalize filename so NF can match output directive
  if   [ -s "${sample_id}.msi" ]; then
    :
  elif [ -s "${sample_id}.msi.txt" ]; then
    mv "${sample_id}.msi.txt" "${sample_id}.msi"
  elif [ -s "${sample_id}" ] && grep -q '^Total_Number_of_Sites' "${sample_id}"; then
    mv "${sample_id}" "${sample_id}.msi"
  elif [ -s "${sample_id}.txt" ] && grep -q '^Total_Number_of_Sites' "${sample_id}.txt"; then
    mv "${sample_id}.txt" "${sample_id}.msi"
  else
    echo "ERROR: Could not find msisensor-pro summary; got:" >&2
    ls -l >&2
    echo "==== ${sample_id}.log (tail) ====" >&2
    tail -n 200 "${sample_id}.log" >&2 || true
    exit 1
  fi

  # Final sanity for Nextflow
  test -s "${sample_id}.msi"
  """
}
