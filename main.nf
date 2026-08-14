#!/usr/bin/env nextflow

// Modules to include
include { CLEANINPUTS               } from './modules/cleaninputs.nf'
include { QUALITYCONTROL            } from './modules/qualitycontrol.nf'
include { CREATE_BREAKS             } from './modules/create_breaks.nf'
include { DESEQ_BASIC               } from './modules/deseq_basic.nf'
include { EDGER_BASIC               } from './modules/edger_basic.nf'
include { WILCOXON_BASIC            } from './modules/wilcoxon_basic.nf'
include { ENET_BASIC                } from './modules/enet_basic.nf'
include { SVC_BASIC                 } from './modules/svc_basic.nf'
//include { SVM_BASIC } from './modules/svm_basic.nf'
include { DATA_PERMUTE              } from './modules/data_permute.nf'
include { DESEQ_PERMUTE             } from './modules/deseq_permute.nf'
include { EDGER_PERMUTE             } from './modules/edger_permute.nf'
include { WILCOXON_PERMUTE          } from './modules/wilcoxon_permute.nf'
include { ENET_PERMUTE              } from './modules/enet_permute.nf'
include { SVC_PERMUTE               } from './modules/svc_permute.nf'
//include { SVM_PERMUTE } from './modules/svm_permute.nf'
include { PERMUTE_PLOTS             } from './modules/permute_plots.nf'
include { PERMUTE_HISTS             } from './modules/permute_hists.nf'
include { INTER_OBSERVER_BASIC      } from './modules/inter_observer_basic.nf'
include { DESEQ_QUASI               } from './modules/deseq_quasi.nf'
include { DESEQ_DATA_QUASI          } from './modules/deseq_data_quasi.nf'
include { EDGER_QUASI               } from './modules/edger_quasi.nf'
include { EDGER_DATA_QUASI          } from './modules/edger_data_quasi.nf'
include { WILCOXON_QUASI            } from './modules/wilcoxon_quasi.nf'
include { WILCOXON_DATA_QUASI       } from './modules/wilcoxon_data_quasi.nf'
//include { SVM_DATA_QUASI } from './modules/svm_data_quasi.nf'
//include { SVM_QUASI } from './modules/svm_quasi.nf'
include { QUASI_PLOTS               } from './modules/quasi_plots.nf'
include { DATA_FULLSYNTH            } from './modules/data_fullsynth.nf'
include { DESEQ_FULLSYNTH           } from './modules/deseq_fullsynth.nf'
include { EDGER_FULLSYNTH           } from './modules/edger_fullsynth.nf'
include { WILCOXON_FULLSYNTH        } from './modules/wilcoxon_fullsynth.nf'
include { ENET_FULLSYNTH            } from './modules/enet_fullsynth.nf'
include { SVC_FULLSYNTH             } from './modules/svc_fullsynth.nf'
include { SVC_FULLSYNTH_POSTPROCESS } from './modules/svc_fullsynth_postprocess.nf'
include { FULLSYNTH_PLOTS           } from './modules/fullsynth_plots.nf'
include { RESAMPLE_ANALYSIS         } from './modules/resample_analysis.nf'
include { RESAMPLE_SUMMARISE        } from './modules/resample_summarise.nf'

//exit 1, 'DEBUG'

// Set number of permutations for fakey datasets
params.nperms    = 7

// By default, run the synthetic data step
params.synthstep = true

// If this is a debug run, pare down the data passed to the ML methods because these take a long time to run
params.pareml    = false
params.mlstep    = true

// Keep the provisional elastic-net branch independently opt-in
params.enetstep  = false

//def cutline = params.samplenum.intdiv(3)
//def breaks = [cutline,cutline*2,params.samplenum]

workflow {

  // Validate inputs
  if (params.samplesheet) {
    ch_samplesheet = file(params.samplesheet)
  }
  else {
    exit(1, 'Input samplesheet not specified!')
  }
  if (params.countsframe) {
    ch_countsframe = file(params.countsframe)
  }
  else {
    exit(1, 'Input count frame not specified!')
  }
  if (params.reflevel) {
    ch_reflevel = params.reflevel
  }
  else {
    exit(1, 'Reference level not specified!')
  }

  println("Sample sheet: " + ch_samplesheet)
  println("Counts matrix: " + ch_countsframe)
  println("Reference level: " + ch_reflevel)

  // Set a vector of depths across which to generate fullsynth datasets
  depths = channel.from(1E+5, 1E+6, 1E+7)
  // Provide a placeholder path when ML outputs are disabled
  dummypath = channel.fromPath("0")

  //Initial data QC and cleanup
  CLEANINPUTS(ch_samplesheet, ch_countsframe, ch_reflevel)
  // Reuse the SVC normalization sidecars for baseline and permuted fits
  svc_library_sizes = CLEANINPUTS.out.svc_library_sizes
  svc_library_genes = CLEANINPUTS.out.svc_library_genes
  // Reuse the pre-filter library sizes for observed elastic-net datasets
  enet_library_sizes = CLEANINPUTS.out.svc_library_sizes
  QUALITYCONTROL(CLEANINPUTS.out.cleaned)
  CREATE_BREAKS(ch_samplesheet)
  // Emit each generated reference/alternative subsample pair
  breaks = CREATE_BREAKS.out.splitCsv(header: false)


  //Basic analysis for each method
  DESEQ_BASIC(CLEANINPUTS.out.cleaned)
  EDGER_BASIC(CLEANINPUTS.out.cleaned)
  WILCOXON_BASIC(CLEANINPUTS.out.cleaned)
  if (params.enetstep) {
    ENET_BASIC(CLEANINPUTS.out.cleaned, enet_library_sizes)
  }
  if (params.mlstep) {
    SVC_BASIC(
      CLEANINPUTS.out.cleaned,
      svc_library_sizes,
      svc_library_genes,
    )
  }

  ////Permutation analysis with no true signal
  perms = channel.from(1..params.nperms)
  DATA_PERMUTE(CLEANINPUTS.out.cleaned, perms)
  // For DESeq
  DESEQ_PERMUTE(DATA_PERMUTE.out.frames)
  // For edgeR
  EDGER_PERMUTE(DATA_PERMUTE.out.frames)
  // For Wilcoxon
  WILCOXON_PERMUTE(DATA_PERMUTE.out.frames)
  if (params.enetstep) {
    // For stability-selected elastic net
    ENET_PERMUTE(DATA_PERMUTE.out.frames, enet_library_sizes)
  }
  if (params.mlstep) {
    // For SVC
    SVC_PERMUTE(
      DATA_PERMUTE.out.frames,
      svc_library_sizes,
      svc_library_genes,
    )
  }
  // Combine outputs and plot
  PERMUTE_PLOTS(
    DESEQ_BASIC.out.table,
    DESEQ_PERMUTE.out.outfile.collect(),
    EDGER_BASIC.out,
    EDGER_PERMUTE.out.outfile.collect(),
    WILCOXON_BASIC.out,
    WILCOXON_PERMUTE.out.outfile.collect(),
    params.mlstep ? SVC_BASIC.out.table : dummypath,
    params.mlstep ? SVC_PERMUTE.out.nFeatures.collect() : dummypath,
    params.enetstep ? ENET_BASIC.out.summary : dummypath,
    params.enetstep
      ? ENET_PERMUTE.out.summary.map { _replicate, summary -> summary }.collect()
      : dummypath,
  )
  // Compare baseline DEG calls across methods
  INTER_OBSERVER_BASIC(
    DESEQ_BASIC.out.table,
    DESEQ_BASIC.out.dds,
    EDGER_BASIC.out,
    WILCOXON_BASIC.out,
    params.mlstep
      ? SVC_BASIC.out.table
      : dummypath,
    params.enetstep
      ? ENET_BASIC.out.table
      : dummypath,
  )

  // Optionally assess end-to-end robustness to fewer biological replicates
  // or lower count depth without crossing the two perturbation grids.
  resampling_configurations = []
  if (params.sample_subsampling || params.depth_subsampling) {
    repeats_text = params.subsample_repeats.toString()
    if (!(repeats_text ==~ /[1-9][0-9]*/)) {
      exit(1, 'subsample_repeats must be a positive integer')
    }
    subsample_repeats_count = repeats_text.toInteger()

    seed_text = params.subsample_seed.toString()
    if (!(seed_text ==~ /[0-9]+/) || seed_text.toBigInteger() > 2147483646) {
      exit(1, 'subsample_seed must be an integer between 0 and 2147483646')
    }
    subsample_seed_value = seed_text.toLong()

    if (params.sample_subsampling) {
      sample_target_tokens = params.subsample_sizes.toString()
        .split(',')*.trim().findAll { token -> token }
      if (!sample_target_tokens ||
          sample_target_tokens.any { token -> !(token ==~ /[1-9][0-9]*/) }) {
        exit(1, 'subsample_sizes must be a comma-separated list of positive integers')
      }
      if (sample_target_tokens.any { token -> token.toBigInteger() > 2147483646 }) {
        exit(1, 'subsample_sizes contains a value above 2147483646')
      }
      sample_targets = sample_target_tokens
        .collect { token -> token.toInteger() }.unique().sort()
      sample_targets.each { target ->
        (1..subsample_repeats_count).each { replicate ->
          def generation_seed = ((subsample_seed_value + target.toLong() * 1009L +
            replicate.toLong() * 9176L) % 2147483646L).toInteger() + 1
          resampling_configurations.add(tuple(
            'sample', replicate, target, generation_seed
          ))
        }
      }
    }

    if (params.depth_subsampling) {
      depth_target_tokens = params.subsample_depths.toString()
        .split(',')*.trim().findAll { token -> token }
      if (!depth_target_tokens ||
          depth_target_tokens.any { token -> !(token ==~ /[1-9][0-9]*/) }) {
        exit(1, 'subsample_depths must be a comma-separated list of positive integers')
      }
      if (depth_target_tokens.any { token -> token.toBigInteger() > 2147483646 }) {
        exit(1, 'subsample_depths contains a value above 2147483646')
      }
      depth_targets = depth_target_tokens
        .collect { token -> token.toInteger() }.unique().sort()
      depth_targets.each { target ->
        (1..subsample_repeats_count).each { replicate ->
          def generation_seed = ((subsample_seed_value + 1000003L +
            target.toLong() * 1009L + replicate.toLong() * 9176L) %
            2147483646L).toInteger() + 1
          resampling_configurations.add(tuple(
            'depth', replicate, target, generation_seed
          ))
        }
      }
    }

    resampling_tasks = channel.fromList(resampling_configurations)
    RESAMPLE_ANALYSIS(
      resampling_tasks,
      ch_samplesheet,
      ch_countsframe,
      ch_reflevel,
    )
    RESAMPLE_SUMMARISE(
      RESAMPLE_ANALYSIS.out.summaries.collect(),
      RESAMPLE_ANALYSIS.out.genes.collect(),
      RESAMPLE_ANALYSIS.out.allocations.collect(),
      RESAMPLE_ANALYSIS.out.libraries.collect(),
      INTER_OBSERVER_BASIC.out.membership,
    )
  }

  // Summarise per-gene selection across null permutations
  PERMUTE_HISTS(
    DESEQ_PERMUTE.out.outfile.collect(),
    EDGER_PERMUTE.out.outfile.collect(),
    WILCOXON_PERMUTE.out.outfile.collect(),
    params.mlstep ? SVC_PERMUTE.out.outfile.collect() : dummypath,
    params.enetstep
      ? ENET_PERMUTE.out.stable.map { _replicate, stable -> stable }.collect()
      : dummypath,
    INTER_OBSERVER_BASIC.out.deseq_poorfits,
  )


  // Quasi-permutation analysis with partial true signal retained
  DESEQ_DATA_QUASI(CLEANINPUTS.out.cleaned, perms.combine(breaks))
  DESEQ_QUASI(DESEQ_DATA_QUASI.out)
  EDGER_DATA_QUASI(CLEANINPUTS.out.cleaned, perms.combine(breaks))
  EDGER_QUASI(EDGER_DATA_QUASI.out)
  WILCOXON_DATA_QUASI(CLEANINPUTS.out.cleaned, perms.combine(breaks))
  WILCOXON_QUASI(WILCOXON_DATA_QUASI.out)
  if (params.mlstep) {
  }
  // Combine outputs and plot
  QUASI_PLOTS(
    DESEQ_QUASI.out.collect(),
    EDGER_QUASI.out.collect(),
    WILCOXON_QUASI.out.collect(),
    dummypath,
  )

  // Allow skipping of this step
  if (params.synthstep) {
    // Synthetic analysis, allowing us to define known true DEGs
    DATA_FULLSYNTH(perms.combine(breaks).combine(depths), ch_countsframe)

    DESEQ_FULLSYNTH(DATA_FULLSYNTH.out)
    EDGER_FULLSYNTH(DATA_FULLSYNTH.out)
    WILCOXON_FULLSYNTH(DATA_FULLSYNTH.out)
    if (params.enetstep) {
      ENET_FULLSYNTH(DATA_FULLSYNTH.out)
    }
    if (params.mlstep) {
      SVC_FULLSYNTH(DATA_FULLSYNTH.out)
      SVC_FULLSYNTH_POSTPROCESS(SVC_FULLSYNTH.out.result)
    }

    FULLSYNTH_PLOTS(
      DESEQ_FULLSYNTH.out.collect(),
      EDGER_FULLSYNTH.out.collect(),
      WILCOXON_FULLSYNTH.out.collect(),
      params.mlstep ? SVC_FULLSYNTH_POSTPROCESS.out.collect() : dummypath,
      params.enetstep ? ENET_FULLSYNTH.out.result.collect() : dummypath,
    )
  }
}
