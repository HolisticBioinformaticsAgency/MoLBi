process CNVKIT_REF {
  tag { ref_fa.baseName }
  container 'quay.io/biocontainers/cnvkit:0.9.12--pyhdfd78af_1'
  env = [ 'MPLCONFIGDIR': '/tmp' ]

  publishDir { "${params.outdir_abs}/reference/cnvkit_ref" }, mode: 'copy'

  input:
  path ref_fa
  path bed

  output:
  path "flat_ref.cnn", emit: refcnn
  // also expose the staged FASTA for downstream use
  path ref_fa,         emit: fasta

  script:
  """
  set -euo pipefail
  cnvkit.py reference \
    -f ${ref_fa} \
    -t ${bed} \
    -o flat_ref.cnn
  """
}
