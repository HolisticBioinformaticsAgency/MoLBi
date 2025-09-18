process FASTQC {
  tag "${sample_id}"
  container = 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
  publishDir "${params.outdir}/fastqc", mode: 'copy'

  input:
  tuple val(sample_id), path(reads)

  output:
  path "*_fastqc.html", emit: html
  path "*_fastqc.zip",  emit: zip

  script:
  """
  fastqc -q -o ./ ${reads.join(" ")}
  """
}

