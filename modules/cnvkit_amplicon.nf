process CNVKIT_AMPLICON {
  tag "${subject}_${sample_id}"
  container 'quay.io/biocontainers/cnvkit:0.9.12--pyhdfd78af_1'
  env = [ 'MPLCONFIGDIR': '/tmp' ]

  publishDir { "${params.outdir_abs}/${subject}/cnvkit" }, mode: 'copy'

  input:
  tuple val(subject), val(sample_id), path(bam), path(bai), path(ref_fa), path(refcnn)

  output:
  tuple val(subject),
        val(sample_id),
        path("${sample_id}.cnr"),
        path("${sample_id}.cns"),
        path("${sample_id}.cnv.vcf"),
        emit: cnv

  script:
  """
  set -euo pipefail

  cnvkit.py batch ${bam} \
    -r ${refcnn} \
    -m amplicon \
    -p ${task.cpus} \
    -d .

  CNR=\$(ls *.cnr | head -n1)
  CNS=\$(ls *.cns | head -n1)

  mv "\$CNR" ${sample_id}.cnr
  mv "\$CNS" ${sample_id}.cns

  # Export CNV calls as VCF (may be empty in some panels)
  cnvkit.py export vcf ${sample_id}.cns -f ${ref_fa} -o ${sample_id}.cnv.vcf || touch ${sample_id}.cnv.vcf
  """
}
