process VEP_ANNOTATE {
  tag "${subject}:${id}"
  container 'quay.io/biocontainers/ensembl-vep:115--pl5321h2a3209d_0'

  // per-subject outputs
  publishDir { "${params.outdir_abs}/${subject}/vep" }, mode: 'copy'

  input:
  // from wiring: (subject, id, vcf, ref_fa)
  tuple val(subject), val(id), path(vcf), path(ref_fa)

  output:
  tuple val(subject), val(id), path("${id}_vep.vcf"),       emit: vepvcf
  tuple val(subject), val(id), path("${id}_vep_stats.txt"), emit: stats

  script:
  """
  set -euo pipefail
  vep \
    --offline --cache \
    --species ${params.vep_species} \
    --assembly ${params.vep_assembly} \
    --dir_cache ${params.vep_cache_dir} \
    --fasta ${ref_fa} \
    --vcf --force_overwrite --format vcf \
    --buffer_size 50 \
    --symbol --numbers --mane --hgvs --total_length \
    --canonical --biotype --protein --ccds --domains --uniprot \
    --af --af_gnomad --max_af --variant_class \
    --pick \
    --input_file ${vcf} \
    --output_file ${id}_vep.vcf \
    --stats_text --stats_file ${id}_vep_stats.txt
  """
}
