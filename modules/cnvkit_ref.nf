process CNVKIT_REF {
  tag { ref_fa.baseName }
  container 'quay.io/biocontainers/cnvkit:0.9.12--pyhdfd78af_1'
  env = [ 'MPLCONFIGDIR': '/tmp' ]

  // Global reference outputs (not per-subject)
  publishDir { "${params.outdir_abs}/reference/cnvkit_ref" }, mode: 'copy'

  input:
  path ref_fa
  path bed

  output:
  path "flat_ref_cnn", emit: refcnn
  // expose the staged FASTA for downstream CNVkit steps
  path ref_fa,         emit: fasta

  script:
  """
  cnvkit.py reference \
    -f ${ref_fa} \
    -t ${bed} \
    -o flat_ref_cnn
  """
}
