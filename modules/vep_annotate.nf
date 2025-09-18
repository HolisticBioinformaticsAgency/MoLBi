process VEP_ANNOTATE {
  tag "${sample_id}"
  container = 'quay.io/biocontainers/ensembl-vep:115--pl5321h2a3209d_0'
  publishDir "${params.outdir}/vep", mode: 'copy'

  input:
  tuple val(sample_id), path(vcf), path(ref_fa)

  output:
  tuple val(sample_id), path("${sample_id}.vep.vcf"), emit: vepvcf
  path "${sample_id}.vep.stats.txt",                  emit: stats

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
    --output_file ${sample_id}.vep.vcf \
    --stats_text --stats_file ${sample_id}.vep.stats.txt
  """
}
