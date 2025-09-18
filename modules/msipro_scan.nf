process MSIPRO_SCAN {
  tag { ref_fa.baseName }
  container = 'quay.io/biocontainers/msisensor-pro:1.3.0--hd979922_1'
  publishDir "${params.outdir}/msisensor", mode: 'copy'

  input:
  path ref_fa

  output:
  path "microsat.sites.tsv", emit: sites
  path ref_fa,               emit: fasta

  script:
  """
  msisensor-pro scan \
    -d ${ref_fa} \
    -o microsat.sites.tsv \
    -p ${task.cpus}
  """
}
