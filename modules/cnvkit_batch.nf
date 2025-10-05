process CNVKIT_BATCH {
  tag "${subject}_${case_id}"
  container 'quay.io/biocontainers/cnvkit:0.9.12--pyhdfd78af_1'
  env = [
    'MPLCONFIGDIR': '/tmp',
    'MPLBACKEND': 'Agg',
    'XDG_CACHE_HOME': '/tmp',
    'HOME': '/tmp'
  ]

  publishDir { "${pub_base}/cnvkit" }, mode: 'copy'

  input:
  tuple val(pub_base), val(subject), val(case_id), val(mode),
        val(main_id),
        path(main_bam, stageAs: 'case.bam'),
        path(main_bai, stageAs: 'case.bam.bai'),
        val(ctrl_id),
        path(ctrl_bam, stageAs: 'control.bam'),
        path(ctrl_bai, stageAs: 'control.bam.bai'),
        path(ref_fa), path(bed), path(refcnn)

  output:
  tuple val(pub_base), val(subject), val(case_id),
        path("${case_id}.cnr"),
        path("${case_id}.cns"),
        emit: cnv
  path "${case_id}-scatter.pdf", optional: true, emit: scatter
  path "${case_id}-diagram.pdf", optional: true, emit: diagram

  script:
  """
  set -euo pipefail

  if [ "${mode}" = "paired" ]; then
    # Build reference on-the-fly from matched normal (no -r)
    cnvkit.py batch case.bam \\
      --normal control.bam \\
      --targets ${bed} \\
      --fasta ${ref_fa} \\
      -p ${task.cpus} \\
      -d .
  else
    # Single-sample: use prebuilt CNN (no --fasta/--targets)
    cnvkit.py batch case.bam \\
      -r ${refcnn} \\
      -p ${task.cpus} \\
      -d .
  fi

  CNR=\$(ls *.cnr | head -n1)
  CNS=\$(ls *.cns | head -n1)
  mv "\$CNR" ${case_id}.cnr
  mv "\$CNS" ${case_id}.cns

  cnvkit.py scatter ${case_id}.cnr -s ${case_id}.cns --segment-color none -g '' -o ${case_id}-scatter.pdf || true
  cnvkit.py diagram ${case_id}.cnr -s ${case_id}.cns -o ${case_id}-diagram.pdf || true
  """
}
