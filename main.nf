nextflow.enable.dsl=2

// ---------- Params  ----------
if( !params.containsKey('reads')        ) params.reads        = "data/*_{1,2}.fastq.gz"
if( !params.containsKey('reference')    ) params.reference    = null
if( !params.containsKey('bed')          ) params.bed          = null
if( !params.containsKey('outdir')       ) params.outdir       = "results"
if( !params.containsKey('min_af')       ) params.min_af       = 0.05
if( !params.containsKey('vardict_mode') ) params.vardict_mode = 'single'  // 'single' | 'paired'

// SnpEff
if( !params.containsKey('snpeff_genome') ) params.snpeff_genome = "GRCh38.99"
if( !params.containsKey('snpeff_data')   ) params.snpeff_data   = ""

// VEP
if( !params.containsKey('vep_cache_dir') ) params.vep_cache_dir = ""
if( !params.containsKey('vep_assembly')  ) params.vep_assembly  = "GRCh38"
if( !params.containsKey('vep_species')   ) params.vep_species   = "homo_sapiens"

// pyTMB annotation choice
if( !params.containsKey('pytmb_annot') ) params.pytmb_annot = "vep"  // "vep" or "snpeff"

// Defaults
if( !params.containsKey('multiqc_extra_args') ) params.multiqc_extra_args = ''

// Compute an absolute outdir once, reuse everywhere
params.outdir_abs = params.outdir.startsWith('/') \
  ? params.outdir \
  : "${projectDir}/${params.outdir}"

// ---------- Sanity ----------
if( !params.reference ) exit 1, "ERROR: --reference FASTA is required"
if( !params.bed )       exit 1, "ERROR: --bed regions BED is required"

// ---------- Includes ----------
include { FASTQC            } from './modules/fastqc.nf'
include { BWA_INDEX         } from './modules/bwa_index.nf'
include { REF_FAIDX         } from './modules/ref_faidx.nf'
include { ALIGN_SAM         } from './modules/align_sam.nf'
include { SORT_INDEX        } from './modules/sort_index.nf'
include { ADD_OR_REPLACE_READ_GROUPS } from './modules/add_or_replace_read_groups.nf'
include { DEDUP_MARKDUPS    } from './modules/dedup_markdups.nf'

include { HSMETRICS } from './modules/hsmetrics.nf'

include { POLYSOLVER } from './modules/polysolver.nf'

include { VARDICT_SINGLE_RAW } from './modules/vardict_single_raw.nf'
include { VARDICT_PAIRED_RAW } from './modules/vardict_paired_raw.nf'
include { VARDICT_TO_VCF     } from './modules/vardict_to_vcf.nf'

include { CNVKIT_REF       } from './modules/cnvkit_ref.nf'
include { CNVKIT_AMPLICON  } from './modules/cnvkit_amplicon.nf'
include { CNVKIT_GENES     } from './modules/cnvkit_genes.nf'

include { MSIPRO_SCAN      } from './modules/msipro_scan.nf'
include { MSIPRO_MSI       } from './modules/msipro_msi.nf'   

include { VEP_ANNOTATE     } from './modules/vep_annotate.nf'
include { SNPEFF_ANNOTATE  } from './modules/snpeff_annotate.nf'
include { CLINVAR_ANNOTATE } from './modules/clinvar_annotate.nf'

include { SOMATIC_FILTER } from './modules/somatic_filter.nf'

include { PYTMB } from './modules/pytmb.nf'

include { MULTIQC } from './modules/multiqc.nf'

// ---------- Channels ----------
Channel.fromFilePairs(params.reads).set { ch_reads }
Channel.of( file(params.reference) ).set { ch_ref_fa }
Channel.of( file(params.bed)       ).set { ch_bed }

// ---------- Workflow ----------
workflow {

  // QC
  FASTQC( ch_reads )

  // Build index & faidx
  BWA_INDEX( ch_ref_fa )
  ch_idxdir  = BWA_INDEX.out.idxdir
  ch_ref_pas = BWA_INDEX.out.fasta

  REF_FAIDX( ch_ref_pas )
  ch_ref_fai = REF_FAIDX.out.fai
  ch_ref_fa2 = REF_FAIDX.out.fasta

  // Align
  ch_align_in = ch_reads.combine(ch_idxdir).combine(ch_ref_pas)
  ch_sam = ALIGN_SAM( ch_align_in )

  // Sort & index (existing)
  ch_bam_sorted = SORT_INDEX( ch_sam )

  // Ensure read groups exist
  ch_bam_rg = ADD_OR_REPLACE_READ_GROUPS( ch_bam_sorted )

  // Remove duplicates on RG-fixed BAMs
  (ch_bam, ch_dedup_metrics) = DEDUP_MARKDUPS( ch_bam_rg )

  // HsMetrics 
  HSMETRICS( ch_bam, ch_bed, ch_ref_fa2, ch_ref_fai )

  // HLA typing
  POLYSOLVER( ch_bam )

  // Tumour–Normal pairing (expects IDs ending with N/T)
  def ch_pairs = ch_bam
    .map { id, bam, bai ->
      def m = (id =~ /(.+?)([NT])$/)
      if( !m ) throw new IllegalArgumentException("Sample id '${id}' must end with 'N' or 'T'")
      tuple(m[0][1], tuple(m[0][2], bam, bai))  // (base, (type, bam, bai))
    }
    .groupTuple()
    .map { base, recs ->
      def t = recs.find { it[0] == 'T' }
      def n = recs.find { it[0] == 'N' }
      if( !t || !n ) { log.warn "Skipping case ${base}: missing ${t ? 'N' : 'T'} mate"; return null }
      tuple(base, t[1], t[2], n[1], n[2])
    }
    .filter { it != null }

  // VarDict inputs
  def ch_vardict_in_single = ch_bam.combine(ch_ref_fa2).combine(ch_ref_fai).combine(ch_bed)
  def ch_vardict_in_paired = ch_pairs.combine(ch_ref_fa2).combine(ch_ref_fai).combine(ch_bed)

  def ch_variants_vcf
  if( params.vardict_mode == 'paired' ) {
    log.info "VarDict mode: paired (tumour–normal)"
    def ch_tsv = VARDICT_PAIRED_RAW( ch_vardict_in_paired ).map{ id, tsv -> tuple(id, 'paired', tsv) }
    ch_variants_vcf = VARDICT_TO_VCF( ch_tsv )
  } else {
    log.info "VarDict mode: single-sample"
    def ch_tsv = VARDICT_SINGLE_RAW( ch_vardict_in_single ).map{ id, tsv -> tuple(id, 'single', tsv) }
    ch_variants_vcf = VARDICT_TO_VCF( ch_tsv )
  }

  // CNVkit reference
  (ch_cnvref_cnn, ch_cnvref_fa) = CNVKIT_REF( ch_ref_fa2, ch_bed )

  // CNVkit per sample
  ch_cnvkit_in = ch_bam.combine(ch_cnvref_fa).combine(ch_cnvref_cnn)
  ch_cnvkit = CNVKIT_AMPLICON( ch_cnvkit_in )

  // CNV per-gene
  CNVKIT_GENES( ch_cnvkit.map{ id, cnr, cns, vcf -> tuple(id, cnr, cns) } )

  // MSI
  (ch_msisites, ch_msiref_fa) = MSIPRO_SCAN( ch_ref_fa2 )
  def ch_msipro_in = ch_pairs.combine(ch_msisites).combine(ch_ref_fa2)
  MSIPRO_MSI( ch_msipro_in )

  // --- Annotation ---

  // ch_variants_vcf emits: tuple(id, vcf)

  // Build (id, vcf, ref_fa) using the param FASTA as a staged file
  def ch_vep_in = ch_variants_vcf.map { sid, vcf ->
    tuple(sid, vcf, file(params.reference))
  }

  // Call the process with the CHANNEL
  def ch_vep    = VEP_ANNOTATE( ch_vep_in )
  
  // Build the single input tuple for SnpEff: (id, vcf, mane_dir, clinvar_vcf)
  def ch_snpeff_in = ch_variants_vcf.map { sid, vcf ->
    tuple(
      sid,
      vcf,
      file("${params.refs_dir}/GRCh38.mane.1.2.refseq"),
      file("${params.refs_dir}/clinvar_20250601.vcf.gz")
    )
  }

  // singletons for MANE dir and ClinVar file (next to main.nf under ./refs)
  def ch_mane_dir    = Channel.of( file("${params.refs_dir}/GRCh38.mane.1.2.refseq") )
  def ch_clinvar_vcf = Channel.of( file("${params.refs_dir}/clinvar_20250601.vcf.gz") )

  // 1) SnpEff (MANE) → emits (id, *.snpeff.core.vcf)
  def ch_snpeff_core_in = ch_variants_vcf
    .combine(ch_mane_dir)                      // (id, vcf) + mane
  def ch_snpeff_core = SNPEFF_ANNOTATE( ch_snpeff_core_in )

  // 2) ClinVar with SnpSift → emits (id, *.snpeff.vcf)
  def ch_clinvar_in = ch_snpeff_core
    .combine(ch_clinvar_vcf)                   // (id, core.vcf) + clinvar
  def ch_snpeff = CLINVAR_ANNOTATE( ch_clinvar_in )

  // --- Pick one annotated stream based on params.pytmb_annot ---
  def ch_annot =
    (params.pytmb_annot == 'snpeff') ? ch_snpeff :
    (params.pytmb_annot == 'vep')    ? ch_vep    :
                                       ch_vep    // default

  // --- Somatic-only filter (run once on the chosen annotation) ---
  def ch_somatic = SOMATIC_FILTER( ch_annot )         // emits: tuple(id, somatic.vcf)

  // --- TMB on the filtered VCFs ---
  PYTMB(
    ch_somatic,
    file(params.pytmb_db_config),
    file(params.pytmb_var_config)
  )

      // ---- MultiQC barrier & run ----
  // HsMetrics files (flatten tuple to path)
  ch_hsmetrics_files = HSMETRICS.out.hs.map { sid, hsfile -> hsfile }

  // Merge with MarkDuplicates metrics
  ch_multiqc_inputs = ch_dedup_metrics.mix(ch_hsmetrics_files).collect()

  // MULTIQC
  MULTIQC( ch_multiqc_inputs )

}