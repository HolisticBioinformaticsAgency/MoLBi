process VEP_ANNOTATE {

  tag "${subject}:${id}"
  container 'quay.io/biocontainers/ensembl-vep:115--pl5321h2a3209d_0'

  input:
  // (pub_base, subject, id, vcf, ref_fa)
  tuple val(pub_base), val(subject), val(id), path(vcf), path(ref_fa)

  output:
  // Main VEP VCF + stats
  tuple val(pub_base), val(subject), val(id),
        path("${id}_vardict_VEP_${params.vep_species}_${params.vep_assembly}_v${params.vep_version}.vcf"),
        emit: vepvcf

  tuple val(pub_base), val(subject), val(id),
        path("${id}_vardict_VEP_${params.vep_species}_${params.vep_assembly}_v${params.vep_version}_stats.txt"),
        emit: stats

  // ProteinSeqs plugin outputs — publish to a sibling dir ProSeqs/
  publishDir { "${pub_base}/ProSeqs" }, mode: 'copy', pattern: "*.fa"

  tuple val(pub_base), val(subject), val(id), path("${id}.ref.fa"), emit: prot_ref
  tuple val(pub_base), val(subject), val(id), path("${id}.mut.fa"), emit: prot_mut

  script:
  """
  set -euo pipefail

  OUT_VCF="${id}_vardict_VEP_${params.vep_species}_${params.vep_assembly}_v${params.vep_version}.vcf"
  OUT_STATS="${id}_vardict_VEP_${params.vep_species}_${params.vep_assembly}_v${params.vep_version}_stats.txt"

  vep \\
    --offline --cache \\
    --species ${params.vep_species} \\
    --assembly ${params.vep_assembly} \\
    --dir_cache ${params.vep_cache_dir} \\
    --fasta ${ref_fa} \\
    --vcf --force_overwrite --format vcf \\
    --buffer_size 50 \\
    --symbol --numbers --mane --hgvs --total_length \\
    --canonical --biotype --protein --ccds --domains --uniprot \\
    --af --af_gnomad --max_af --variant_class \\
    --pick \\
    --plugin ProteinSeqs,${id}.ref.fa,${id}.mut.fa \\
    --input_file ${vcf} \\
    --output_file \${OUT_VCF} \\
    --stats_text --stats_file \${OUT_STATS}
  """
}
