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

  // ---- Param check ---------
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

  // ---------- Build per-subject groups ----------
  def by_subject = ch_bam_meta
    .map { sub, sample, status, sex, bam, bai -> tuple(sub, tuple(sample, status, bam, bai)) }
    .groupTuple()

  // ---------- Compute cases ----------
  def ch_cases = by_subject.flatMap { sub, recs ->
    def tumors  = recs.findAll { it[1] == 'tumor'  }
                   .collect { [ it[0], it[2], it[3], 'tumor'  ] }
    def normals = recs.findAll { it[1] == 'normal' }
                   .collect { [ it[0], it[2], it[3], 'normal' ] }

    def out = []
    if( tumors && normals ) {
      tumors.each { t ->
        normals.each { n ->
          def case_id = "${t[0]}_${n[0]}"
          def samples = [ t, n ]
          out << tuple(sub, case_id, 'paired', samples)
        }
      }
    } else {
      recs.each { r ->
        def sample_id = r[0]; def bam = r[2]; def bai = r[3]; def status = r[1]
        def samples = [ [ sample_id, bam, bai, status ] ]
        out << tuple(sub, sample_id, 'single', samples)
      }
    }
    out
  }

  // ---------- Decide publishing base ----------
  def ch_case_counts = ch_cases
    .map { sub, case_id, mode, samples -> tuple(sub, 1) }
    .groupTuple()
    .map { sub, ones -> tuple(sub, ones.size()) }

  def ch_cases_pub = ch_cases
    .join(ch_case_counts)
    .map { sub, case_id, mode, samples, ncases ->
      def base = (ncases > 1) ? "${params.outdir_abs}/${sub}/${case_id}" : "${params.outdir_abs}/${sub}"
      tuple(sub, case_id, mode, samples, base)
    }

  // ---------- POLYSOLVER ----------
  POLYSOLVER(
    ch_cases_pub.flatMap { sub, case_id, mode, samples, pub_base ->
      samples.collect { s ->
        tuple(pub_base, sub, case_id, s[0], s[1], s[2])
      }
    }
  )

  // ---------- VarDict ----------
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

  def ch_vardict_single_vcf = VARDICT_SINGLE( ch_vardict_single_in )
    .map { pub_base, sub, id, vcf -> tuple(pub_base, sub, id, 'single', vcf) }

  def ch_vardict_paired_vcf = VARDICT_PAIRED( ch_vardict_paired_in )
    .map { pub_base, sub, case_id, vcf -> tuple(pub_base, sub, case_id, 'paired', vcf) }

  def ch_variants_vcf = ch_vardict_single_vcf.mix( ch_vardict_paired_vcf )

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

  def ch_cnv, ch_scatter, ch_diagram
  (ch_cnv, ch_scatter, ch_diagram) = CNVKIT_BATCH( ch_cnv_jobs )

  CNVKIT_GENES(
    ch_cnv.map { pub_base, sub, case_id, cnr, cns -> tuple(pub_base, sub, case_id, cnr, cns) }
  )

  // ---------- MSI ----------
  (ch_msisites, ch_msiref_fa) = MSIPRO_SCAN( ch_ref_fa )

  // Paired MSI (msisensor-pro msi)
  def ch_msi_jobs_paired = ch_cases_pub
    .filter { sub, case_id, mode, samples, pub_base -> mode == 'paired' }
    .combine( ch_msisites )
    .combine( ch_msiref_fa )
    .flatMap { sub, case_id, mode, samples, pub_base, sites_tsv, ref_fa ->
      def tumors  = samples.findAll { it[3] == 'tumor' }
      def normals = samples.findAll { it[3] == 'normal' }
      tumors.collectMany { t ->
        normals.collect { n ->
          // (pub_base, subject, case_id, t_id, t_bam, t_bai, n_id, n_bam, n_bai, sites, fasta)
          tuple(pub_base, sub, case_id, t[0], t[1], t[2], n[0], n[1], n[2], sites_tsv, ref_fa)
        }
      }
    }

  // Single tumour MSI (msisensor-pro pro)
  def ch_msi_jobs_single_tumor = ch_cases_pub
    .filter { sub, case_id, mode, samples, pub_base ->
      mode == 'single' && samples[0][3] == 'tumor'
    }
    .combine( ch_msisites )
    .combine( ch_msiref_fa )
    .map { sub, case_id, mode, samples, pub_base, sites_tsv, ref_fa ->
      def t = samples[0]
      // (pub_base, subject, case_id, t_id, t_bam, t_bai, sites, fasta)
      tuple(pub_base, sub, case_id, t[0], t[1], t[2], sites_tsv, ref_fa)
    }

  MSIPRO_PAIRED ( ch_msi_jobs_paired )
  MSIPRO_SINGLE ( ch_msi_jobs_single_tumor )

  // ---------- Annotation ----------
  def ch_vep_in = ch_variants_vcf.combine(ch_ref_fa)
    .map { pub_base, sub, id, mode, vcf, ref_fa ->
      tuple(pub_base, sub, id, file(vcf), file(ref_fa))
    }
  def ch_vep_vcf, ch_vep_stats
  (ch_vep_vcf, ch_vep_stats) = VEP_ANNOTATE( ch_vep_in )

  def ch_mane_dir    = Channel.of( file(params.mane_dir) )
  def ch_snpeff_core = SNPEFF_ANNOTATE(
    ch_variants_vcf.combine(ch_mane_dir).map { pub_base, sub, id, mode, vcf, mane ->
      tuple(pub_base, sub, id, mode, vcf, mane)
    }
  )
  def ch_clinvar_vcf = Channel.of( file(params.clinvar_vcf) )
  def ch_clinvar_out = CLINVAR_ANNOTATE(
    ch_snpeff_core.combine(ch_clinvar_vcf).map { pub_base, sub, id, mode, core_vcf, clin ->
      tuple(pub_base, sub, id, mode, core_vcf, clin)
    }
  )

  def ch_annot_for_all = (params.pytmb_annot == 'snpeff')
    ? ch_clinvar_out.map { pub_base, sub, id, mode, vcf -> tuple(pub_base, sub, id, mode, vcf) }
    : ch_vep_vcf.combine( ch_variants_vcf.map{ pb, s, i, m, v -> tuple(s,i,m) } )
                 .map { vep_pb_sub_id_vcf, sub_id_mode ->
                   def (pb, s1, i1, vepvcf) = vep_pb_sub_id_vcf
                   def (s2, i2, mode)      = sub_id_mode
                   assert s1==s2 && i1==i2
                   tuple(pb, s1, i1, mode, vepvcf)
                 }

  // Stop normal-only after annotation
  def ch_cases_has_tumor = ch_cases_pub.map { sub, case_id, mode, samples, pub_base ->
    tuple("${sub}::${case_id}", sub, case_id, samples.any{ it[3]=='tumor' })
  }
  def ch_annot_keyed    = ch_annot_for_all.map { pub, sub, id, mode, vcf -> tuple("${sub}::${id}", pub, sub, id, mode, vcf) }
  def ch_anno_with_flag = ch_annot_keyed.join( ch_cases_has_tumor )
    .map { key, pub, sub, id, mode, vcf, sub2, case_id, hasTumor ->
      tuple(pub, sub, id, mode, vcf, hasTumor)
    }

  def ch_for_tumor = ch_anno_with_flag
    .filter { pub, sub, id, mode, vcf, hasTumor -> hasTumor }
    .map    { pub, sub, id, mode, vcf, hasTumor -> tuple(pub, sub, id, mode, vcf) }

  def ch_somatic_filtered_vcf = SOMATIC_FILTER( ch_for_tumor )

  def ch_cases_keyed = ch_cases_pub.map { sub, case_id, mode, samples, pub_base ->
    tuple("${sub}::${case_id}", sub, case_id, mode, samples, pub_base)
  }
  def ch_somatic_key = ch_somatic_filtered_vcf.map { pub, sub, id, vcf ->
    tuple("${sub}::${id}", pub, sub, id, vcf)
  }

  def ch_tmb_joined = ch_somatic_key.join( ch_cases_keyed )
    .map { key, pub, sub, id, vcf, sub2, case_id, mode, samples, pub_base ->
      def tumorRec   = samples.find { it[3] == 'tumor' }
      def tumorLabel = tumorRec ? tumorRec[0] : null
      tuple(pub_base ?: pub, sub, case_id, tumorLabel, vcf)
    }

  def ch_tmb_jobs = ch_tmb_joined.filter { pub, sub, case_id, tumorLabel, vcf -> tumorLabel != null }

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

workflow.onComplete {
    def outdir = params.outdir_abs
    def cmd = """
      set -euo pipefail
      shopt -s nullglob

      for d in "${outdir}"/* ; do
        [ -d "\$d" ] || continue
        base="\$(basename "\$d")"
        # Skip non-subject folders at top level
        case "\$base" in
          reference|multiqc) continue ;;
        esac

        dest="\${d%/}/vcf"
        mkdir -p "\$dest"

        # Gather all VCFs & indexes anywhere under the subject dir,
        # but don't recurse into the destination vcf/ we're filling.
        find "\$d" -type f \\
          \\( -name "*.vcf" -o -name "*.vcf.gz" -o -name "*.vcf.tbi" -o -name "*.vcf.gz.tbi" -o -name "*.csi" \\) \\
          -not -path "\$dest/*" \\
          -exec cp -f {} "\$dest/" \\;
      done
    """
    def p = ["bash","-lc", cmd].execute()
    p.consumeProcessOutput(System.out, System.err)
    def rc = p.waitFor()
    if( rc != 0 ) log.warn "VCF gather post-step exited with code ${rc}"
}
