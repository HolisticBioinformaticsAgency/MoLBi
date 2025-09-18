process CNVKIT_AMPLICON {
  tag "${sample_id}"
  container 'quay.io/biocontainers/cnvkit:0.9.12--pyhdfd78af_1'
  env = [ 'MPLCONFIGDIR': '/tmp' ]

  publishDir "${params.outdir}/cnvkit", mode: 'copy'

  input:
  tuple val(sample_id), path(bam), path(bai), path(ref_fa), path(refcnn)

  output:
  tuple val(sample_id),
        path("${sample_id}.cnr"),
        path("${sample_id}.cns"),
        path("${sample_id}.cnv.vcf"), emit: cnv

  script:
  """
  cnvkit.py batch ${bam} \
    -r ${refcnn} \
    -m amplicon \
    -p ${task.cpus} \
    -d .

  CNR=\$(ls *.cnr | head -n1); CNS=\$(ls *.cns | head -n1)
  mv "\$CNR" ${sample_id}.cnr
  mv "\$CNS" ${sample_id}.cns

  cnvkit.py export vcf ${sample_id}.cns -f ${ref_fa} -o ${sample_id}.cnv.vcf || touch ${sample_id}.cnv.vcf
  """
}
