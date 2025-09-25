process CNVKIT_REF {
  tag { ref_fa.baseName }
  container 'quay.io/biocontainers/cnvkit:0.9.12--pyhdfd78af_1'
  env = [ 'MPLCONFIGDIR': '/tmp' ]

<<<<<<< HEAD
  // Global reference outputs (not per-subject)
=======
>>>>>>> f12105e (Pipeline moved to vh83, dropped filtering flags for sort_index)
  publishDir { "${params.outdir_abs}/reference/cnvkit_ref" }, mode: 'copy'

  input:
  path ref_fa
  path bed

  output:
<<<<<<< HEAD
  path "flat_ref_cnn", emit: refcnn
  // expose the staged FASTA for downstream CNVkit steps
=======
  path "flat_ref.cnn", emit: refcnn
  // also expose the staged FASTA for downstream use
>>>>>>> f12105e (Pipeline moved to vh83, dropped filtering flags for sort_index)
  path ref_fa,         emit: fasta

  script:
  """
<<<<<<< HEAD
  cnvkit.py reference \
    -f ${ref_fa} \
    -t ${bed} \
    -o flat_ref_cnn
=======
  set -euo pipefail
  cnvkit.py reference \
    -f ${ref_fa} \
    -t ${bed} \
    -o flat_ref.cnn
>>>>>>> f12105e (Pipeline moved to vh83, dropped filtering flags for sort_index)
  """
}
