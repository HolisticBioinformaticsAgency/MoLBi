// modules/init_params.nf
nextflow.enable.dsl=2

/**
 * INIT_PARAMS
 * - Sets defaults ONLY when a key is missing (so YAML/CLI override them)
 * - Derives params.outdir_abs as an absolute path for stable publishDir targets
 */
workflow INIT_PARAMS {

  // ---- Defaults (only if user didn't provide them via -params-file or CLI) ----
  if( !params.containsKey('samplesheet')        ) params.samplesheet         = "./refs/samplesheet.csv"
  if( !params.containsKey('reference')          ) params.reference           = null
  if( !params.containsKey('reference_fai')      ) params.reference_fai       = null
  if( !params.containsKey('bed')                ) params.bed                 = null
  if( !params.containsKey('outdir')             ) params.outdir              = "results"
  if( !params.containsKey('min_af')             ) params.min_af              = 0.05
  if( !params.containsKey('vardict_mode')       ) params.vardict_mode        = 'paired'   // 'paired' | 'single'
  if( !params.containsKey('multiqc_extra_args') ) params.multiqc_extra_args  = ''

  // Optional tool params (silence warnings; modules may still use their own defaults)
  if( !params.containsKey('polysolver_race')       ) params.polysolver_race       = 'Caucasian'
  if( !params.containsKey('polysolver_build')      ) params.polysolver_build      = 'hg38'
  if( !params.containsKey('polysolver_emit_vcf')   ) params.polysolver_emit_vcf   = 0
  if( !params.containsKey('polysolver_fastqtype')  ) params.polysolver_fastqtype  = 'STDFQ'
  if( !params.containsKey('polysolver_insertcalc') ) params.polysolver_insertcalc = 0

  // ---- Derived values (need Nextflow context, not YAML) ----
  // Absolute outdir (stable for publishDir closures)
  params.outdir_abs = params.outdir?.startsWith('/') \
    ? params.outdir \
    : "${projectDir}/${params.outdir}"

  // (Optional) You can add early sanity checks here if you want, but keeping
  // this module strictly for defaults + derivations keeps it lightweight.
}

