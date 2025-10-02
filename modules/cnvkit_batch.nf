process CNVKIT_BATCH {
  tag "${subject}_${tumor_id}"
  container 'quay.io/biocontainers/cnvkit:0.9.12--pyhdfd78af_1'
  env = [
    'MPLCONFIGDIR': '/tmp',
    'MPLBACKEND': 'Agg',
    'XDG_CACHE_HOME': '/tmp',
    'HOME': '/tmp'
  ]

  publishDir { "${params.outdir_abs}/${subject}/cnvkit" }, mode: 'copy'

  input:
  // (subject, tumor_id, t_bam, t_bai, normal_id, n_bam, n_bai, ref_fa, bed)
  tuple val(subject),
        val(tumor_id),  path(t_bam), path(t_bai),
        val(normal_id), path(n_bam), path(n_bai),
        path(ref_fa),    path(bed)

  output:
  tuple val(subject),
        val(tumor_id),
        path("${tumor_id}.cnr"),
        path("${tumor_id}.cns"),
        emit: cnv
  path "${tumor_id}-scatter.pdf",  optional: true, emit: scatter
  path "${tumor_id}-diagram.pdf",  optional: true, emit: diagram

  script:
  """
  set -euo pipefail

  # CNVkit tumor–normal: build reference from matched normal (no -r)
  cnvkit.py batch ${t_bam} \\
    --normal ${n_bam} \\
    --targets ${bed} \\
    --fasta ${ref_fa} \\
    -p ${task.cpus} \\
    -d .

  # Standardize filenames
  CNR=\$(ls *.cnr | head -n1)
  CNS=\$(ls *.cns | head -n1)
  mv "\$CNR" ${tumor_id}.cnr
  mv "\$CNS" ${tumor_id}.cns

  # Make plots from CNR/CNS (don’t fail pipeline if plotting hiccups)
  cnvkit.py scatter ${tumor_id}.cnr -s ${tumor_id}.cns -o ${tumor_id}-scatter.pdf  || true
  cnvkit.py diagram ${tumor_id}.cnr -s ${tumor_id}.cns -o ${tumor_id}-diagram.pdf || true
  """
}
