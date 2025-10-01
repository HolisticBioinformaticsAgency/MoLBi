nextflow.enable.dsl=2

// Absolute outdir (stable for publishDir closures)
params.outdir_abs = params.outdir?.startsWith('/') \
  ? params.outdir \
  : "${projectDir}/${params.outdir}"

// -------------------- Includes --------------------
include { INIT_PARAMS }                    from './modules/init_params.nf'
include { FASTQC }                         from './modules/fastqc.nf'
include { ALIGN_AND_SORT }                 from './modules/align_and_sort.nf'
include { DEDUP_MARKDUPS }                 from './modules/dedup_markdups.nf'
include { HSMETRICS }                      from './modules/hsmetrics.nf'
include { POLYSOLVER }                     from './modules/polysolver.nf'
include { VARDICT_SINGLE_RAW }             from './modules/vardict_single_raw.nf'
include { VARDICT_PAIRED_RAW }             from './modules/vardict_paired_raw.nf'
include { VARDICT_TO_VCF }                 from './modules/vardict_to_vcf.nf'
include { CNVKIT_REF }                     from './modules/cnvkit_ref.nf'
include { CNVKIT_AMPLICON }                from './modules/cnvkit_amplicon.nf'
include { CNVKIT_GENES }                   from './modules/cnvkit_genes.nf'
include { MSIPRO_SCAN }                    from './modules/msipro_scan.nf'
include { MSIPRO_MSI }                     from './modules/msipro_msi.nf'
include { VEP_ANNOTATE }                   from './modules/vep_annotate.nf'
include { SNPEFF_ANNOTATE }                from './modules/snpeff_annotate.nf'
include { CLINVAR_ANNOTATE }               from './modules/clinvar_annotate.nf'
include { SOMATIC_FILTER }                 from './modules/somatic_filter.nf'
include { PYTMB }                          from './modules/pytmb.nf'
include { MULTIQC }                        from './modules/multiqc.nf'

// -------------------- Samplesheet → Channels --------------------
def REQUIRED = ['sample','subject','status','fastq_1','fastq_2']

Channel
  .fromPath(params.samplesheet)
  .splitCsv(header:true)
  .map { row ->
    def missing = REQUIRED.findAll { k ->
      !row.containsKey(k) || row[k] == null || row[k].toString().trim() == ''
    }
    if( missing ) {
      throw new IllegalArgumentException("Samplesheet row has missing/empty field(s): ${missing} → ${row}")
    }
    tuple(
      row.sample as String,
      row.subject as String,
      (row.status as String).toLowerCase(),   // 'tumor' | 'normal'
      file(row.fastq_1),
      file(row.fastq_2),
      ((row.sex ?: 'NA') as String)
    )
  }
  .set { ch_sheet }  // (sample, subject, status, R1, R2, sex)

// Reads for QC/BWA: (subject, sample, [R1,R2])
ch_reads = ch_sheet.map { sample, subject, status, r1, r2, sex ->
  tuple(subject, sample, [ r1, r2 ])
}

// Sample metadata for joins: (sample, subject, status, sex)
ch_meta  = ch_sheet.map { sample, subject, status, r1, r2, sex ->
  tuple(sample, subject, status, sex)
}

// Reference & BED as singletons
Channel.of( file(params.reference)     ).set { ch_ref_fa }
Channel.of( file(params.reference_fai) ).set { ch_ref_fai }
Channel.of( file(params.bed)           ).set { ch_bed }

// Also pass the ORIGINAL absolute path string for the reference (so we can find sidecars)
Channel
  .of( params.reference )
  .map { new File(it as String).getAbsolutePath() }
  .set { ch_ref_src_abs }

// -------------------- Workflow --------------------
workflow {

  // ---- Param check ---------
  INIT_PARAMS()

  // ---------- QC ----------
  FASTQC( ch_reads )

  // ---------- Align + Sort (no SAM on disk) ----------
  ch_align_in = ch_reads
    .combine(ch_ref_fa)
    .combine(ch_ref_src_abs)

  // Emits (subject, sample_id, *.hq.sorted.bam, *.hq.sorted.bam.bai)
  ch_bam_sorted = ALIGN_AND_SORT( ch_align_in )

  // ---------- Dedup ----------
  (ch_bam, ch_dedup_metrics) = DEDUP_MARKDUPS( ch_bam_sorted )

  // ---------- Attach metadata (subject/status/sex) ----------
  ch_bam_meta = ch_bam
    .map  { sub, sample, bam, bai -> tuple(sample, tuple(sub, bam, bai)) }
    .join ( ch_meta )
    .map  { sample, left, subject, status, sex ->
      tuple(subject, sample, status, sex, left[1], left[2])
    }

  // ---------- HsMetrics ----------
  def ch_hs_in = ch_bam
    .combine(ch_bed)
    .combine(ch_ref_fa)
    .combine(ch_ref_fai)
    .map { sub, sample, bam, bai, bed, ref_fa, ref_fai ->
      tuple(sub, sample, bam, bai, bed, ref_fa, ref_fai)
    }

  HSMETRICS( ch_hs_in )

  // ---------- HLA typing (POLYSOLVER) ----------
  def ch_bam_for_hla_in =
  ch_bam_meta
    .map { sub, sample, status, sex, bam, bai -> tuple(sub, sample, bam, bai) }

  POLYSOLVER( ch_bam_for_hla_in )


  // ---------- Pair tumour/normal by subject ----------
  def ch_pairs = ch_bam_meta
    .map { sub, sample, status, sex, bam, bai -> tuple(sub, tuple(status, sample, bam, bai)) }
    .groupTuple()
    .map { sub, recs ->
      def t = recs.find { it[0] == 'tumor' }
      def n = recs.find { it[0] == 'normal' }
      if( !t || !n ) { log.warn "Skipping ${sub}: missing ${t ? 'normal' : 'tumor'}"; return null }
      tuple(sub, t[1], t[2], t[3], n[1], n[2], n[3])
    }
    .filter { it != null }

  // ---------- VarDict ----------
  def ch_vardict_in_single = ch_bam_meta
    .combine(ch_ref_fa)
    .combine(ch_ref_fai)
    .combine(ch_bed)
    .map { sub, sample, status, sex, bam, bai, ref_fa, ref_fai, bed ->
      tuple(sub, sample, bam, bai, ref_fa, ref_fai, bed)
    }

  def ch_vardict_in_paired = ch_pairs
    .combine(ch_ref_fa)
    .combine(ch_ref_fai)
    .combine(ch_bed)
    .map { sub, t_sample, t_bam, t_bai, n_sample, n_bam, n_bai, ref_fa, ref_fai, bed ->
      tuple(sub, t_sample, t_bam, t_bai, n_sample, n_bam, n_bai, ref_fa, ref_fai, bed)
    }

  def ch_vcf_in
  if( params.vardict_mode == 'paired' ) {
    log.info "VarDict mode: paired (tumour–normal)"
    ch_vcf_in = VARDICT_PAIRED_RAW( ch_vardict_in_paired )
      .map { sampleId, tsv -> tuple(sampleId, 'paired', null, null, null, tsv, sampleId) }
  } else {
    log.info "VarDict mode: single-sample"
    ch_vcf_in = VARDICT_SINGLE_RAW( ch_vardict_in_single )
      .map { sampleId, tsv -> tuple(sampleId, 'single', null, null, sampleId, tsv, sampleId) }
  }
  def ch_variants_vcf = VARDICT_TO_VCF( ch_vcf_in )

  // ---------- CNVkit ----------
  (ch_cnvref_cnn, ch_cnvref_fa) = CNVKIT_REF( ch_ref_fa, ch_bed )
  def ch_cnvkit_in = ch_bam_meta
    .combine(ch_cnvref_fa)
    .combine(ch_cnvref_cnn)
    .map { sub, sample, status, sex, bam, bai, ref_fa, cnn ->
      tuple(sub, sample, bam, bai, ref_fa, cnn)
    }
  def ch_cnvkit = CNVKIT_AMPLICON( ch_cnvkit_in )
  CNVKIT_GENES( ch_cnvkit.map { sub, sample, cnr, cns, vcf -> tuple(sub, sample, cnr, cns) } )

  // ---------- MSI ----------
  (ch_msisites, ch_msiref_fa) = MSIPRO_SCAN( ch_ref_fa )
  def ch_msipro_in = ch_pairs
    .combine( ch_msisites )
    .combine( ch_msiref_fa )
    .map { sub, t_sample, t_bam, t_bai, n_sample, n_bam, n_bai, sites_tsv, ref_fa ->
      tuple(sub, t_sample, t_bam, t_bai, n_sample, n_bam, n_bai, sites_tsv, ref_fa)
    }
  MSIPRO_MSI( ch_msipro_in )

  // ---------- Annotation ----------
  def ch_vep_in = ch_variants_vcf.map { sub, sid, vcf -> tuple(sub, sid, vcf, file(params.reference)) }
  def ch_vep    = VEP_ANNOTATE( ch_vep_in )

  def ch_mane_dir    = Channel.of( file(params.mane_dir) )
  def ch_clinvar_vcf = Channel.of( file(params.clinvar_vcf) )

  def ch_mane_annot = SNPEFF_ANNOTATE(
    ch_variants_vcf.combine(ch_mane_dir)
                   .map { sub, sid, vcf, mane -> tuple(sub, sid, vcf, mane) }
  )
  def ch_snpeff = CLINVAR_ANNOTATE(
    ch_mane_annot.combine(ch_clinvar_vcf)
                  .map { sub, sid, core_vcf, clin -> tuple(sub, sid, core_vcf, clin) }
  )

  def ch_annot = (params.pytmb_annot == 'snpeff') ? ch_snpeff :
                 (params.pytmb_annot == 'vep')    ? ch_vep    :
                                                    ch_vep

  def ch_somatic_filtered_vcf = SOMATIC_FILTER( ch_annot )

  // ----------- TMB -----------------
  PYTMB(
    ch_somatic_filtered_vcf.map { sub, sid, vcf -> tuple(sub, sid, file(vcf)) },
    file(params.pytmb_db_config),
    file(params.pytmb_var_config)
  )

  // ---------- MultiQC ----------
  def ch_hsmetrics_files = HSMETRICS.out.hs.map { sub, sid, hs -> hs }
  def ch_multiqc_inputs  = ch_dedup_metrics.mix(ch_hsmetrics_files).collect()
  MULTIQC( ch_multiqc_inputs )
}

workflow.onComplete {
    def outdir = params.outdir_abs ?: "${projectDir}/results"
    def cmd = """
      set -euo pipefail
      shopt -s nullglob

      for d in "${outdir}"/* ; do
        [ -d "\$d" ] || continue
        dest="\${d%/}/vcf"
        mkdir -p "\$dest"

        # Copy VCFs and indexes from immediate subdirs, but skip the vcf/ folder itself
        find "\$d" -maxdepth 2 -type f \\
          \\( -name "*.vcf.gz" -o -name "*.vcf.tbi" -o -name "*.vcf.gz.tbi" -o -name "*.csi" \\) \\
          -not -path "\$dest/*" \\
          -exec cp -f {} "\$dest/" \\;
      done
    """
    def p = ["bash","-lc", cmd].execute()
    p.consumeProcessOutput(System.out, System.err)
    def rc = p.waitFor()
    if( rc != 0 ) log.warn "VCF gather post-step exited with code ${rc}"
}



