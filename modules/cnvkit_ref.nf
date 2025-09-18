process CNVKIT_REF {
  tag { ref_fa.baseName }
  container 'quay.io/biocontainers/cnvkit:0.9.12--pyhdfd78af_1'
  env = [ 'MPLCONFIGDIR': '/tmp' ]

  publishDir "${params.outdir}/cnvkit_ref", mode: 'copy'

  input:
  path ref_fa
  path bed

  output:
  path "flat_ref.cnn", emit: refcnn
  path ref_fa,         emit: fasta

  script:
  """
  cnvkit.py reference \
    -f ${ref_fa} \
    -t ${bed} \
    -o flat_ref.cnn
  """
}
