#!/usr/bin/env nextflow

// Modules to include
include { CLEANINPUTS } from './modules/cleaninputs.nf'
include { QUALITYCONTROL } from './modules/qualitycontrol.nf'
include { DESEQ_BASIC } from './modules/deseq_basic.nf'
include { EDGER_BASIC } from './modules/edger_basic.nf'
include { WILCOXON_BASIC } from './modules/wilcoxon_basic.nf'
include { SVC_BASIC } from './modules/svc_basic.nf'
include { DATA_PERMUTE } from './modules/data_permute.nf'
include { DESEQ_PERMUTE } from './modules/deseq_permute.nf'
include { EDGER_PERMUTE } from './modules/edger_permute.nf'
include { WILCOXON_PERMUTE } from './modules/wilcoxon_permute.nf'
include { SVC_PERMUTE } from './modules/svc_permute.nf'
include { PERMUTE_PLOTS } from './modules/permute_plots.nf'
include { PERMUTE_HISTS } from './modules/permute_hists.nf'
include { DESEQ_QUASI } from './modules/deseq_quasi.nf'
include { DESEQ_DATA_QUASI } from './modules/deseq_data_quasi.nf'
include { EDGER_QUASI } from './modules/edger_quasi.nf'
include { EDGER_DATA_QUASI } from './modules/edger_data_quasi.nf'
include { WILCOXON_QUASI } from './modules/wilcoxon_quasi.nf'
include { WILCOXON_DATA_QUASI } from './modules/wilcoxon_data_quasi.nf'
include { DATA_FULLSYNTH } from './modules/data_fullsynth.nf'
include { DESEQ_FULLSYNTH } from './modules/deseq_fullsynth.nf'
include { EDGER_FULLSYNTH } from './modules/edger_fullsynth.nf'
include { WILCOXON_FULLSYNTH } from './modules/wilcoxon_fullsynth.nf'
include { SVC_FULLSYNTH } from './modules/svc_fullsynth.nf'
include { SVC_FULLSYNTH_POSTPROCESS } from './modules/svc_fullsynth_postprocess.nf'
include { FULLSYNTH_PLOTS } from './modules/fullsynth_plots.nf'

//exit 1, 'DEBUG'

// Validate inputs
if (params.samplesheet) { ch_samplesheet = file(params.samplesheet) } else { exit 1, 'Input samplesheet not specified!' }
if (params.countsframe) { ch_countsframe = file(params.countsframe) } else { exit 1, 'Input count frame not specified!' }
if (params.reflevel) { ch_reflevel = params.reflevel } else { exit 1, 'Reference level not specified!' }

// Set number of permutations for fakey datasets
params.nperms = 7

println("Sample sheet: " + ch_samplesheet)
println("Counts matrix: " + ch_countsframe)
println("Reference level: " + ch_reflevel)

//def cutline = params.samplenum.intdiv(3)
//def breaks = [cutline,cutline*2,params.samplenum]

// Set a vector of depths across which to generate fullsynth datasets
depths = Channel.from(1e5, 1e6, 1e7)

process QUASI_PLOTS {

  publishDir "$params.outdir"

  input:
  val deseq_quasis
  val edger_quasis
  val wilcox_quasis

  output:
  path "quasi_plot.pdf"

  """
  #!/usr/bin/env Rscript
  library("tidyverse", quietly = TRUE)
  library("reshape2", quietly = TRUE)

  deseq.inlist = str_remove_all("$deseq_quasis","[\\\\[\\\\] ]") %>% 
    strsplit(split = ",") %>% unlist()

  deseq.quasi.input = 
    sapply(deseq.inlist, read.csv) %>% 
    t() %>% 
    `row.names<-`(NULL) %>% 
    data.frame() %>%
    mutate_all(as.numeric) %>%
    mutate(method = "DESeq2") %>% 
    melt(id.vars = c("method","samplenum"))

  edger.inlist = str_remove_all("$edger_quasis","[\\\\[\\\\] ]") %>% 
    strsplit(split = ",") %>% unlist()

  edger.quasi.input = 
    sapply(edger.inlist, read.csv) %>% 
    t() %>% 
    `row.names<-`(NULL) %>% 
    data.frame() %>%
    mutate_all(as.numeric) %>%
    mutate(method = "edgeR") %>% 
    melt(id.vars = c("method","samplenum"))

  wilcox.inlist = str_remove_all("$wilcox_quasis","[\\\\[\\\\] ]") %>% 
    strsplit(split = ",") %>% unlist()

  wilcox.quasi.input = 
    sapply(wilcox.inlist, read.csv) %>% 
    t() %>% 
    `row.names<-`(NULL) %>% 
    data.frame() %>%
    mutate_all(as.numeric) %>%
    mutate(method = "wilcox") %>% 
    melt(id.vars = c("method","samplenum"))

  gg.quasi.input = rbind(deseq.quasi.input, 
                          edger.quasi.input,
                          wilcox.quasi.input)

  gg.quasi = ggplot(gg.quasi.input, aes(x = method, y = value)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_y_continuous(limits = c(0,1)) +
    labs(x = "Method",y = "Replicates per group") +
    facet_grid(samplenum~variable) +
    theme_bw() +
    theme(strip.background = element_blank())

  ggsave(gg.quasi, 
     filename = "quasi_plot.pdf",
     device = "pdf", bg = "transparent",
     width =  30, height = 20, units = "cm")
  """
}

process CREATE_BREAKS {

  debug true

  input:
  path samplesheet

  output:
  file 'breaks.csv'
  //env refouts, emit: refouts
  //env altouts, emit: altouts

  """
  #!/usr/bin/env Rscript
  library("tidyverse", quietly = TRUE)
  samplesheet = read.csv("$samplesheet")
  reflevel = "$params.reflevel"

  reftot = nrow(subset(samplesheet, phenotype == reflevel))
  alttot = nrow(subset(samplesheet, phenotype != reflevel))

  refbreaks = cut(seq(1,reftot),breaks = 3, right = TRUE) %>% levels() %>%
    str_match("\\\\,[^]]*") %>%
    str_remove(",") %>%
    as.numeric() %>%
    round()

  altbreaks = cut(seq(1,alttot),breaks = 3, right = TRUE) %>% levels() %>%
    str_match("\\\\,[^]]*") %>%
    str_remove(",") %>%
    as.numeric() %>%
    round()

  outbreaks = data.frame(refbreaks, altbreaks)
  write.table(outbreaks, "breaks.csv", sep = ",", row.names = F, col.names = F)
  """
}



workflow {
  //Initial data QC and cleanup
  CLEANINPUTS(ch_samplesheet, ch_countsframe, ch_reflevel)
  QUALITYCONTROL(CLEANINPUTS.out)
  CREATE_BREAKS(ch_samplesheet)
  breaks = CREATE_BREAKS.out.splitCsv(header: false)


  //Basic analysis for each method
  DESEQ_BASIC(CLEANINPUTS.out)
  EDGER_BASIC(CLEANINPUTS.out)
  WILCOXON_BASIC(CLEANINPUTS.out)
  SVC_BASIC(CLEANINPUTS.out)

  ////Permutation analysis with no true signal
  perms = Channel.from(1..params.nperms)
  DATA_PERMUTE(CLEANINPUTS.out, perms)
  // For DESeq
  DESEQ_PERMUTE(DATA_PERMUTE.out.frames)
  // For edgeR
  EDGER_PERMUTE(DATA_PERMUTE.out.frames)
  // For Wilcoxon
  WILCOXON_PERMUTE(DATA_PERMUTE.out.frames)
  // For SVC
  SVC_PERMUTE(DATA_PERMUTE.out)
  // Combine outputs and plot
  PERMUTE_PLOTS(
    DESEQ_BASIC.out.table,
    DESEQ_PERMUTE.out.nDEGs.collect(),
    EDGER_BASIC.out,
    EDGER_PERMUTE.out.nDEGs.collect(),
    WILCOXON_BASIC.out,
    WILCOXON_PERMUTE.out.nDEGs.collect(),
    SVC_BASIC.out.table,
    SVC_PERMUTE.out.nDEGs.collect(),
  )
  PERMUTE_HISTS(
    DESEQ_PERMUTE.out.outfile.collect(),
    EDGER_PERMUTE.out.outfile.collect(),
    WILCOXON_PERMUTE.out.outfile.collect(),
    SVC_PERMUTE.out.outfile.collect()
  )

  // Quasi-permutation analysis with partial true signal retained
  DESEQ_DATA_QUASI(CLEANINPUTS.out, perms.combine(breaks))
  DESEQ_QUASI(DESEQ_DATA_QUASI.out)
  EDGER_DATA_QUASI(CLEANINPUTS.out, perms.combine(breaks))
  EDGER_QUASI(EDGER_DATA_QUASI.out)
  WILCOXON_DATA_QUASI(CLEANINPUTS.out, perms.combine(breaks))
  WILCOXON_QUASI(WILCOXON_DATA_QUASI.out)
  // Combine outputs and plot
  QUASI_PLOTS(
    DESEQ_QUASI.out.collect(),
    EDGER_QUASI.out.collect(),
    WILCOXON_QUASI.out.collect()
  )

  // Synthetic analysis, allowing us to define known true DEGs
  DATA_FULLSYNTH(perms.combine(breaks).combine(depths))

  DESEQ_FULLSYNTH(DATA_FULLSYNTH.out)
  EDGER_FULLSYNTH(DATA_FULLSYNTH.out)
  WILCOXON_FULLSYNTH(DATA_FULLSYNTH.out)
  SVC_FULLSYNTH(DATA_FULLSYNTH.out)
  SVC_FULLSYNTH_POSTPROCESS(SVC_FULLSYNTH.out)

  FULLSYNTH_PLOTS(DESEQ_FULLSYNTH.out.collect(),
                  EDGER_FULLSYNTH.out.collect(),
                  WILCOXON_FULLSYNTH.out.collect(),
                  SVC_FULLSYNTH_POSTPROCESS.out.collect())

}
