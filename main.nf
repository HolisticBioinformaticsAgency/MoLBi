nextflow.enable.dsl=2

// ------------------------------------------------------------
// Absolute outdir (stable for publishDir closures)
params.outdir_abs = params.outdir?.startsWith('/') \
  ? params.outdir \
  : "${projectDir}/${params.outdir ?: 'results'}"

// -------------------- Includes --------------------
include { INIT_PARAMS }                    from './modules/init_params.nf'
include { FASTQC }                         from './modules/fastqc.nf'
include { ALIGN_AND_SORT }                 from './modules/align_and_sort.nf'
include { DEDUP_MARKDUPS }                 from './modules/dedup_markdups.nf'
include { HSMETRICS }                      from './modules/hsmetrics.nf'
include { POLYSOLVER }                     from './modules/polysolver.nf'

include { VARDICT_SINGLE }                 from './modules/vardict.nf'
include { VARDICT_PAIRED }                 from './modules/vardict.nf'

include { CNVKIT_REF }                     from './modules/cnvkit_ref.nf'
include { CNVKIT_BATCH }                   from './modules/cnvkit_batch.nf'
include { CNVKIT_GENES }                   from './modules/cnvkit_genes.nf'

include { MSIPRO_SCAN }                    from './modules/msipro_scan.nf'
include { MSIPRO_PAIRED ; MSIPRO_SINGLE }  from './modules/msipro.nf'

include { VEP_ANNOTATE }                   from './modules/vep_annotate.nf'
include { SNPEFF_ANNOTATE }                from './modules/snpeff_annotate.nf'
include { CLINVAR_ANNOTATE }               from './modules/clinvar_annotate.nf'
include { SOMATIC_FILTER }                 from './modules/somatic_filter.nf'
include { PYTMB }                          from './modules/pytmb.nf'
include { MULTIQC }                        from './modules/multiqc.nf'

// ---- same module, multiple aliases ----
include { BGZIP_TABIX as BGZIP_VARDICT_SINGLE } from './modules/bgzip_tabix.nf'
include { BGZIP_TABIX as BGZIP_VARDICT_PAIRED } from './modules/bgzip_tabix.nf'
include { BGZIP_TABIX as BGZIP_SNPEFF }         from './modules/bgzip_tabix.nf'
include { BGZIP_TABIX as BGZIP_CLINVAR }        from './modules/bgzip_tabix.nf'
include { BGZIP_TABIX as BGZIP_VEP }            from './modules/bgzip_tabix.nf'
include { BGZIP_TABIX as BGZIP_SOMATIC }        from './modules/bgzip_tabix.nf'

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
      (row.status as String).toLowerCase(),
      file(row.fastq_1),
      file(row.fastq_2),
      ((row.sex ?: 'NA') as String)
    )
  }
  .set { ch_sheet }

// Reads for QC/BWA: (subject, sample, [R1,R2])
ch_reads = ch_sheet.map { sample, subject, status, r1, r2, sex ->
  tuple(subject, sample, [ r1, r2 ])
}

// Sample metadata for joins: (sample, subject, status, sex)
ch_meta  = ch_sheet.map { sample, subject, status, r1, r2, sex ->
  tuple(sample, subject, status, sex)
}

// Reference & BED singletons
Channel.of( file(params.reference)     ).set { ch_ref_fa }
Channel.of( file(params.reference_fai) ).set { ch_ref_fai }
Channel.of( file(params.bed)           ).set { ch_bed }

// Also pass original absolute FASTA path (for BWA sidecars)
Channel
  .of( params.reference )
  .map { new File(it as String).getAbsolutePath() }
  .set { ch_ref_src_abs }

// -------------------- Workflow --------------------
workflow {

  INIT_PARAMS()

  // ---------- QC ----------
  FASTQC( ch_reads )

  // ---------- Align + Sort ----------
  ch_align_in   = ch_reads.combine(ch_ref_fa).combine(ch_ref_src_abs)
  ch_bam_sorted = ALIGN_AND_SORT( ch_align_in )

  // ---------- Dedup ----------
  (ch_bam, ch_dedup_metrics) = DEDUP_MARKDUPS( ch_bam_sorted )

  // ---------- Attach metadata ----------
  ch_bam_meta = ch_bam
    .map  { sub, sample, bam, bai -> tuple(sample, tuple(sub, bam, bai)) }
    .join ( ch_meta )
    .map  { sample, left, subject, status, sex ->
      tuple(subject, sample, status, sex, left[1], left[2])
    }

  // ---------- HsMetrics ----------
  def ch_hs_in = ch_bam.combine(ch_bed).combine(ch_ref_fa).combine(ch_ref_fai)
    .map { sub, sample, bam, bai, bed, ref_fa, ref_fai -> tuple(sub, sample, bam, bai, bed, ref_fa, ref_fai) }
  HSMETRICS( ch_hs_in )

  // ---------- Group by subject ----------
  def by_subject = ch_bam_meta
    .map { sub, sample, status, sex, bam, bai -> tuple(sub, tuple(sample, status, bam, bai)) }
    .groupTuple()
    .map { sub, recs -> tuple(sub, recs.sort { a, b -> a[0] <=> b[0] }) }

  // ---------- Compute cases and publishing base ----------
  def ch_cases_pub = by_subject.flatMap { sub, recs ->
    def tumors  = recs.findAll { it[1] == 'tumor'  }.collect { [ it[0], it[2], it[3], 'tumor'  ] }
    def normals = recs.findAll { it[1] == 'normal' }.collect { [ it[0], it[2], it[3], 'normal' ] }

    def ncases = (tumors && normals) ? (tumors.size() * normals.size()) : recs.size()

    def out = []
    if( tumors && normals ) {
      tumors.each { t ->
        normals.each { n ->
          def case_id = "${t[0]}_${n[0]}"
          def base    = (ncases > 1) ? "${params.outdir_abs}/${sub}/${case_id}" : "${params.outdir_abs}/${sub}"
          out << tuple(sub, case_id, 'paired', [t, n], base)
        }
      }
    } else {
      recs.each { r ->
        def case_id = r[0]
        def base    = (ncases > 1) ? "${params.outdir_abs}/${sub}/${case_id}" : "${params.outdir_abs}/${sub}"
        out << tuple(sub, case_id, 'single', [ [ r[0], r[2], r[3], r[1] ] ], base)
      }
    }
    out
  }

  // ---------- POLYSOLVER ----------
  POLYSOLVER(
    ch_cases_pub.flatMap { sub, case_id, mode, samples, pub_base ->
      samples.collect { s ->
        tuple(pub_base, sub, case_id, s[0], s[1], s[2])
      }
    }
  )

  // ---------- VarDict in ----------
  def ch_vardict_single_in = ch_cases_pub
    .filter { sub, case_id, mode, samples, pub_base -> mode == 'single' }
    .combine(ch_ref_fa).combine(ch_ref_fai).combine(ch_bed)
    .map { sub, case_id, mode, samples, pub_base, ref_fa, ref_fai, bed ->
      def s = samples[0]
      tuple(pub_base, sub, s[0], s[1], s[2], ref_fa, ref_fai, bed)
    }

  def ch_vardict_paired_in = ch_cases_pub
    .filter { sub, case_id, mode, samples, pub_base -> mode == 'paired' }
    .combine(ch_ref_fa).combine(ch_ref_fai).combine(ch_bed)
    .flatMap { sub, case_id, mode, samples, pub_base, ref_fa, ref_fai, bed ->
      def tumors  = samples.findAll { it[3] == 'tumor'  }
      def normals = samples.findAll { it[3] == 'normal' }
      tumors.collectMany { t ->
        normals.collect { n ->
          tuple(pub_base, sub, case_id, t[0], t[1], t[2], n[0], n[1], n[2], ref_fa, ref_fai, bed)
        }
      }
    }

  // ---------- VarDict run ----------
  VARDICT_SINGLE( ch_vardict_single_in )
  def ch_vardict_single_vcf_raw = VARDICT_SINGLE.out.vcf
    .map { pub, sub, id, vcf -> tuple(pub, 'vardict', sub, id, 'single', vcf) }

  VARDICT_PAIRED( ch_vardict_paired_in )
  def ch_vardict_paired_vcf_raw = VARDICT_PAIRED.out.vcf
    .map { pub, sub, case_id, vcf -> tuple(pub, 'vardict', sub, case_id, 'paired', vcf) }

  // ---------- bgzip/tabix VarDict ----------
  def ch_vardict_single_vcfgz = BGZIP_VARDICT_SINGLE( ch_vardict_single_vcf_raw ).vcfgz
    .map { pub, bucket, sub, id, mode, vcf -> tuple(pub, sub, id, mode, vcf) }

  def ch_vardict_paired_vcfgz = BGZIP_VARDICT_PAIRED( ch_vardict_paired_vcf_raw ).vcfgz
    .map { pub, bucket, sub, id, mode, vcf -> tuple(pub, sub, id, mode, vcf) }

  def ch_variants_vcf = ch_vardict_single_vcfgz.mix( ch_vardict_paired_vcfgz )

  // ---------- CNVkit ref ----------
  (ch_cnvref_cnn, ch_cnvref_fa) = CNVKIT_REF( ch_ref_fa, ch_bed )

  // ---------- CNVkit per case ----------
  def ch_cnv_jobs = ch_cases_pub
    .combine(ch_cnvref_fa).combine(ch_bed).combine(ch_cnvref_cnn)
    .flatMap { sub, case_id, mode, samples, pub_base, ref_fa, bed, cnn ->
      if( mode == 'paired' ) {
        def tumors  = samples.findAll { it[3] == 'tumor' }
        def normals = samples.findAll { it[3] == 'normal' }
        tumors.collectMany { t ->
          normals.collect { n ->
            tuple(pub_base, sub, case_id, 'paired', t[0], t[1], t[2], n[0], n[1], n[2], ref_fa, bed, cnn)
          }
        }
      } else {
        def s = samples[0]
        [ tuple(pub_base, sub, case_id, 'single',
                s[0], s[1], s[2],
                'NA', file('/dev/null'), file('/dev/null'),
                ref_fa, bed, cnn) ]
      }
    }

  CNVKIT_BATCH( ch_cnv_jobs )
  def ch_cnv     = CNVKIT_BATCH.out.cnv
  def ch_scatter = CNVKIT_BATCH.out.scatter
  def ch_diagram = CNVKIT_BATCH.out.diagram

  CNVKIT_GENES(
    ch_cnv.map { pub_base, sub, case_id, cnr, cns -> tuple(pub_base, sub, case_id, cnr, cns) }
  )

  // ---------- MSI ----------
  MSIPRO_SCAN( ch_ref_fa )
  def ch_msisites  = MSIPRO_SCAN.out.sites
  def ch_msiref_fa = MSIPRO_SCAN.out.fasta

  def ch_msi_jobs_paired = ch_cases_pub
    .filter { sub, case_id, mode, samples, pub_base -> mode == 'paired' }
    .combine( ch_msisites )
    .combine( ch_msiref_fa )
    .flatMap { sub, case_id, mode, samples, pub_base, sites_tsv, ref_fa ->
      def tumors  = samples.findAll { it[3] == 'tumor' }
      def normals = samples.findAll { it[3] == 'normal' }
      tumors.collectMany { t ->
        normals.collect { n ->
          tuple(pub_base, sub, case_id, t[0], t[1], t[2], n[0], n[1], n[2], sites_tsv, ref_fa)
        }
      }
    }

  def ch_msi_jobs_single_tumor = ch_cases_pub
    .filter { sub, case_id, mode, samples, pub_base -> mode == 'single' && samples[0][3] == 'tumor' }
    .combine( ch_msisites )
    .combine( ch_msiref_fa )
    .map { sub, case_id, mode, samples, pub_base, sites_tsv, ref_fa ->
      def t = samples[0]
      tuple(pub_base, sub, case_id, t[0], t[1], t[2], sites_tsv, ref_fa)
    }

  MSIPRO_PAIRED ( ch_msi_jobs_paired )
  MSIPRO_SINGLE ( ch_msi_jobs_single_tumor )

  // ---------- VEP ----------
  def ch_vep_in = ch_variants_vcf.combine(ch_ref_fa)
    .map { pair ->
      def (pub_base, sub, id, mode, vcf, ref_fa) = pair
      tuple(pub_base, sub, id, file(vcf), file(ref_fa))
    }
  VEP_ANNOTATE( ch_vep_in )
  def ch_vep_vcf_raw = VEP_ANNOTATE.out.vepvcf

  def ch_vep_vcf = BGZIP_VEP(
    ch_vep_vcf_raw.map { pb, sub, id, vcf -> tuple(pb, 'vep', sub, id, 'vep', vcf) }
  ).vcfgz.map { pb, bucket, sub, id, mode, vcf -> tuple(pb, sub, id, vcf) }

  def ch_vep_stats = VEP_ANNOTATE.out.stats

  // ---------- SnpEff ----------
  def ch_mane_dir = Channel.of( file(params.mane_dir) )

  SNPEFF_ANNOTATE(
    ch_variants_vcf.combine(ch_mane_dir)
      .map { pair ->
        def (pub_base, sub, id, mode, vcf, mane) = pair
        tuple(pub_base, sub, id, mode, vcf, mane)
      }
  )
  def ch_snpeff_core_raw = SNPEFF_ANNOTATE.out.snpeff

  def ch_snpeff_core = BGZIP_SNPEFF(
    ch_snpeff_core_raw.map { pb, sub, id, mode, vcf -> tuple(pb, 'snpeff_MANE', sub, id, mode, vcf) }
  ).vcfgz.map { pb, bucket, sub, id, mode, vcf -> tuple(pb, sub, id, mode, vcf) }

  // ---------- ClinVar ----------
  def ch_clinvar_vcf = Channel.of( file(params.clinvar_vcf) )

  CLINVAR_ANNOTATE(
    ch_snpeff_core.combine(ch_clinvar_vcf)
      .map { pair ->
        def (pub_base, sub, id, mode, core_vcf, clin) = pair
        tuple(pub_base, sub, id, mode, core_vcf, clin)
      }
  )
  def ch_clinvar_out_raw = CLINVAR_ANNOTATE.out.snpeffvcf

  def ch_clinvar_out = BGZIP_CLINVAR(
    ch_clinvar_out_raw.map { pb, sub, id, mode, vcf -> tuple(pb, 'ClinVar', sub, id, mode, vcf) }
  ).vcfgz.map { pb, bucket, sub, id, mode, vcf -> tuple(pb, sub, id, mode, vcf) }

  // ---------- Choose annotation stream ----------
  def ch_annot_for_all = (params.pytmb_annot == 'snpeff')
    ? ch_clinvar_out.map { pb, sub, id, mode, vcf -> tuple(pb, sub, id, mode, vcf) }
    : ch_vep_vcf    .map { pb, sub, id, vcf ->           tuple(pb, sub, id, 'vep', vcf) }

  // ---------- Mark which cases have tumor ----------
  def ch_cases_has_tumor = ch_cases_pub.map { sub, case_id, mode, samples, pub_base ->
    tuple("${sub}::${case_id}", sub, case_id, samples.any { it[3] == 'tumor' })
  }

  def ch_annot_keyed = ch_annot_for_all.map { pb, sub, id, mode, vcf ->
    tuple("${sub}::${id}", pb, sub, id, mode, vcf)
  }

  def ch_anno_with_flag = ch_annot_keyed.join( ch_cases_has_tumor )
    .map { pair ->
      def (key, pb, sub, id, mode, vcf, sub2, case_id, hasTumor) = pair
      tuple(pb, sub, id, mode, vcf, hasTumor)
    }

  def ch_for_tumor = ch_anno_with_flag
    .filter { pb, sub, id, mode, vcf, hasTumor -> hasTumor }
    .map    { pb, sub, id, mode, vcf, hasTumor -> tuple(pb, sub, id, mode, vcf) }

  // ---------- Somatic filter ----------
  SOMATIC_FILTER( ch_for_tumor )
  def ch_somatic_filtered_vcf_raw = SOMATIC_FILTER.out.somatic

  def ch_somatic_filtered_vcf = BGZIP_SOMATIC(
    ch_somatic_filtered_vcf_raw.map { pb, sub, id, vcf -> tuple(pb, 'somatic', sub, id, 'somatic', vcf) }
  ).vcfgz.map { pb, bucket, sub, id, mode, vcf -> tuple(pb, sub, id, vcf) }

  // ---------- TMB ----------
  def ch_cases_keyed = ch_cases_pub.map { sub, case_id, mode, samples, pub_base ->
    tuple("${sub}::${case_id}", sub, case_id, mode, samples, pub_base)
  }

  def ch_somatic_key = ch_somatic_filtered_vcf.map { pb, sub, id, vcf ->
    tuple("${sub}::${id}", pb, sub, id, vcf)
  }

  def ch_tmb_joined = ch_somatic_key.join( ch_cases_keyed )
    .map { pair ->
      def (key, pb, sub, id, vcf, sub2, case_id, mode, samples, pub_base) = pair
      def tumorRec   = samples.find { it[3] == 'tumor' }
      def tumorLabel = tumorRec ? tumorRec[0] : null
      tuple(pub_base ?: pb, sub, case_id, tumorLabel, vcf)
    }

  def ch_tmb_jobs = ch_tmb_joined.filter { pb, sub, case_id, tumorLabel, vcf -> tumorLabel != null }

  PYTMB(
    ch_tmb_jobs,
    file(params.pytmb_db_config),
    file(params.pytmb_var_config)
  )

  // ---------- MultiQC ----------
  def ch_hsmetrics_files = HSMETRICS.out.hs.map { sub, sid, hs -> hs }
  def ch_multiqc_inputs  = ch_dedup_metrics.mix(ch_hsmetrics_files).collect()
  MULTIQC( ch_multiqc_inputs )
}

