#!/usr/bin/env nextflow

// Modules to include
include { CLEANINPUTS               } from './modules/cleaninputs.nf'
include { QUALITYCONTROL            } from './modules/qualitycontrol.nf'
include { CREATE_BREAKS             } from './modules/create_breaks.nf'
include { DESEQ_BASIC               } from './modules/deseq_basic.nf'
include { EDGER_BASIC               } from './modules/edger_basic.nf'
include { WILCOXON_BASIC            } from './modules/wilcoxon_basic.nf'
include { SVC_BASIC                 } from './modules/svc_basic.nf'
//include { SVM_BASIC } from './modules/svm_basic.nf'
include { DATA_PERMUTE              } from './modules/data_permute.nf'
include { DESEQ_PERMUTE             } from './modules/deseq_permute.nf'
include { EDGER_PERMUTE             } from './modules/edger_permute.nf'
include { WILCOXON_PERMUTE          } from './modules/wilcoxon_permute.nf'
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
include { SVC_FULLSYNTH             } from './modules/svc_fullsynth.nf'
include { SVC_FULLSYNTH_POSTPROCESS } from './modules/svc_fullsynth_postprocess.nf'
include { FULLSYNTH_PLOTS           } from './modules/fullsynth_plots.nf'

//exit 1, 'DEBUG'

// Set number of permutations for fakey datasets
params.nperms    = 7

// By default, run the synthetic data step
params.synthstep = true

// If this is a debug run, pare down the data passed to the ML methods because these take a long time to run
params.pareml    = false
params.mlstep    = true

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
  depths = Channel.from(1E+5, 1E+6, 1E+7)
  dummypath = Channel.fromPath("0")

  //Initial data QC and cleanup
  CLEANINPUTS(ch_samplesheet, ch_countsframe, ch_reflevel)
  QUALITYCONTROL(CLEANINPUTS.out)
  CREATE_BREAKS(ch_samplesheet)
  breaks = CREATE_BREAKS.out.splitCsv(header: false)


  //Basic analysis for each method
  DESEQ_BASIC(CLEANINPUTS.out)
  EDGER_BASIC(CLEANINPUTS.out)
  WILCOXON_BASIC(CLEANINPUTS.out)
  if (params.mlstep) {
    SVC_BASIC(CLEANINPUTS.out)
  }

  ////Permutation analysis with no true signal
  perms = Channel.from(1..params.nperms)
  DATA_PERMUTE(CLEANINPUTS.out, perms)
  // For DESeq
  DESEQ_PERMUTE(DATA_PERMUTE.out.frames)
  // For edgeR
  EDGER_PERMUTE(DATA_PERMUTE.out.frames)
  // For Wilcoxon
  WILCOXON_PERMUTE(DATA_PERMUTE.out.frames)
  if (params.mlstep) {
    // For SVC
    SVC_PERMUTE(DATA_PERMUTE.out)
  }
  // Combine outputs and plot
  PERMUTE_PLOTS(
    DESEQ_BASIC.out.table,
    DESEQ_PERMUTE.out.nDEGs.collect(),
    EDGER_BASIC.out,
    EDGER_PERMUTE.out.nDEGs.collect(),
    WILCOXON_BASIC.out,
    WILCOXON_PERMUTE.out.nDEGs.collect(),
    params.mlstep ? SVC_BASIC.out.table : dummypath,
    params.mlstep ? SVC_PERMUTE.out.nDEGs.collect() : dummypath,
  )
  INTER_OBSERVER_BASIC(
    DESEQ_BASIC.out.table,
    DESEQ_BASIC.out.dds,
    EDGER_BASIC.out,
    WILCOXON_BASIC.out,
    params.mlstep
      ? SVC_BASIC.out.table
      : dummypath,
  )
  PERMUTE_HISTS(
    DESEQ_PERMUTE.out.outfile.collect(),
    EDGER_PERMUTE.out.outfile.collect(),
    WILCOXON_PERMUTE.out.outfile.collect(),
    params.mlstep ? SVC_PERMUTE.out.outfile.collect() : dummypath,
    INTER_OBSERVER_BASIC.out.deseq_poorfits,
  )


  // Quasi-permutation analysis with partial true signal retained
  DESEQ_DATA_QUASI(CLEANINPUTS.out, perms.combine(breaks))
  DESEQ_QUASI(DESEQ_DATA_QUASI.out)
  EDGER_DATA_QUASI(CLEANINPUTS.out, perms.combine(breaks))
  EDGER_QUASI(EDGER_DATA_QUASI.out)
  WILCOXON_DATA_QUASI(CLEANINPUTS.out, perms.combine(breaks))
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
    if (params.mlstep) {
      SVC_FULLSYNTH(DATA_FULLSYNTH.out)
      SVC_FULLSYNTH_POSTPROCESS(SVC_FULLSYNTH.out)
    }

    FULLSYNTH_PLOTS(
      DESEQ_FULLSYNTH.out.collect(),
      EDGER_FULLSYNTH.out.collect(),
      WILCOXON_FULLSYNTH.out.collect(),
      params.mlstep ? SVC_FULLSYNTH_POSTPROCESS.out.collect() : dummypath,
    )
  }
}
