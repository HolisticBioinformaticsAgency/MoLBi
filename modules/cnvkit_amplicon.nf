process CNVKIT_AMPLICON {
  tag "${subject}_${sample_id}"
  container 'quay.io/biocontainers/cnvkit:0.9.12--pyhdfd78af_1'
  env = [ 'MPLCONFIGDIR': '/tmp' ]

<<<<<<< HEAD
  // per-subject outputs
=======
>>>>>>> f12105e (Pipeline moved to vh83, dropped filtering flags for sort_index)
  publishDir { "${params.outdir_abs}/${subject}/cnvkit" }, mode: 'copy'

  input:
  tuple val(subject), val(sample_id), path(bam), path(bai), path(ref_fa), path(refcnn)

  output:
  tuple val(subject),
        val(sample_id),
<<<<<<< HEAD
        path("${sample_id}_cnr"),
        path("${sample_id}_cns"),
        path("${sample_id}_cnv.vcf"),
=======
        path("${sample_id}.cnr"),
        path("${sample_id}.cns"),
        path("${sample_id}.cnv.vcf"),
>>>>>>> f12105e (Pipeline moved to vh83, dropped filtering flags for sort_index)
        emit: cnv

  script:
  """
  set -euo pipefail

  cnvkit.py batch ${bam} \
    -r ${refcnn} \
    -m amplicon \
    -p ${task.cpus} \
    -d .

<<<<<<< HEAD
  CNR=\$(ls *.cnr | head -n1); CNS=\$(ls *.cns | head -n1)
  mv "\$CNR" ${sample_id}_cnr
  mv "\$CNS" ${sample_id}_cns

  # Export CNV calls as VCF (may be empty in some panels)
  cnvkit.py export vcf ${sample_id}_cns -f ${ref_fa} -o ${sample_id}_cnv.vcf || touch ${sample_id}_cnv.vcf
=======
  CNR=\$(ls *.cnr | head -n1)
  CNS=\$(ls *.cns | head -n1)

  mv "\$CNR" ${sample_id}.cnr
  mv "\$CNS" ${sample_id}.cns

  # Export CNV calls as VCF (may be empty in some panels)
  cnvkit.py export vcf ${sample_id}.cns -f ${ref_fa} -o ${sample_id}.cnv.vcf || touch ${sample_id}.cnv.vcf
>>>>>>> f12105e (Pipeline moved to vh83, dropped filtering flags for sort_index)
  """
}
