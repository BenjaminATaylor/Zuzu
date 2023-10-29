#!/usr/bin/env nextflow

// Modules to include
include { CLEANINPUTS } from './modules/cleaninputs.nf'
include { QUALITYCONTROL } from './modules/qualitycontrol.nf'
include { DESEQ_BASIC } from './modules/deseq_basic.nf'
include { DATA_PERMUTE } from './modules/data_permute.nf'
include { DESEQ_PERMUTE } from './modules/deseq_permute.nf'
include { PERMUTE_PLOTS } from './modules/permute_plots.nf'
include { PERMUTE_HISTS } from './modules/permute_hists.nf'
include { DESEQ_DATA_QUASI } from './modules/deseq_data_quasi.nf'
include { DESEQ_QUASI } from './modules/deseq_quasi.nf'
include { EDGER_QUASI } from './modules/edger_quasi.nf'
include { EDGER_DATA_QUASI } from './modules/edger_data_quasi.nf'

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

def cutline = params.samplenum.intdiv(3)
def breaks = [cutline,cutline*2,params.samplenum]


process EDGER_BASIC{

  input: 
  tuple path(samplesheet), path(countsframe) 

  output:
  path "edger_table.csv"

  script:
  """
  #!/usr/bin/env Rscript
  library("edgeR") 
  library("tidyverse")

  # Input data
  samplesheet = read.csv("$samplesheet")
  countsframe.clean = read.csv("$countsframe", row.names = 1)

  # Basic edgeR workflow, taken from the documentation
  dds.edge = DGEList(counts=countsframe.clean,group=samplesheet\$phenotype)
  dds.edge = normLibSizes(dds.edge)
  design.edge = model.matrix(~samplesheet\$phenotype)
  dds.edge.est = estimateDisp(dds.edge,design.edge)
  edge.fit = glmQLFit(dds.edge.est,design.edge)
  edge.qlf = glmQLFTest(edge.fit,coef=2)

  # Tabular output (perform FDR correction directly rather than relying on topTags to do it for us)
  edge.out = edge.qlf\$table %>% mutate(padj = p.adjust(PValue, method = "BH")) 
  write.csv(edge.out, row.names = T, file="edger_table.csv")
  """
}

process EDGER_PERMUTE{

  input: 
  tuple path(samplesheet), path(countsframe) 

  output:
  path "edger_table.csv", emit: outfile
  path "nDEGs.RData", emit: nDEGs

  script:
  """
  #!/usr/bin/env Rscript
  library("edgeR") 
  library("tidyverse")

  # Input data
  samplesheet = read.csv("$samplesheet")
  countsframe.clean = read.csv("$countsframe", row.names = 1)

  # Basic edgeR workflow, taken from the documentation
  dds.edge = DGEList(counts=countsframe.clean,group=samplesheet\$phenotype)
  dds.edge = normLibSizes(dds.edge)
  design.edge = model.matrix(~samplesheet\$phenotype)
  dds.edge.est = estimateDisp(dds.edge,design.edge)
  edge.fit = glmQLFit(dds.edge.est,design.edge)
  edge.qlf = glmQLFTest(edge.fit,coef=2)

  # Tabular output (perform FDR correction directly rather than relying on topTags to do it for us)
  edge.out = edge.qlf\$table %>% mutate(padj = p.adjust(PValue, method = "BH")) 

  nDEGs = nrow(subset(edge.out,padj<0.05))
  save(nDEGs,file = "nDEGs.RData")
  write.csv(edge.out, file="edger_table.csv", row.names = T)
  """
}

process QUASI_PLOTS {

  publishDir "$params.outdir"

  input:
  val deseq_quasis
  val edger_quasis

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

  gg.quasi.input = rbind(deseq.quasi.input, edger.quasi.input)

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



workflow {
  //Initial data QC and cleanup
  CLEANINPUTS(ch_samplesheet, ch_countsframe, ch_reflevel)
  QUALITYCONTROL(CLEANINPUTS.out)

  //Basic analysis for each method
  DESEQ_BASIC(CLEANINPUTS.out)
  EDGER_BASIC(CLEANINPUTS.out)

  ////Permutation analysis with no true signal
  perms = Channel.from(1..params.nperms)
  DATA_PERMUTE(CLEANINPUTS.out, perms)
  // For DESeq
  DESEQ_PERMUTE(DATA_PERMUTE.out)
  DESEQ_PERMUTE.out.nDEGs
  // For edgeR
  EDGER_PERMUTE(DATA_PERMUTE.out)
  EDGER_PERMUTE.out.nDEGs
  // Combine outputs and plot
  PERMUTE_PLOTS(
    DESEQ_BASIC.out.table,
    DESEQ_PERMUTE.out.nDEGs.collect(),
    EDGER_BASIC.out,
    EDGER_PERMUTE.out.nDEGs.collect()
  )
  PERMUTE_HISTS(
    DESEQ_PERMUTE.out.outfile.collect(),
    EDGER_PERMUTE.out.outfile.collect()
  )

  //Quasi-permutation analysis with partial true signal retained
  samplenums = Channel.from(breaks)
  //samplenums.subscribe { println "DEBUG: $it" }

  DESEQ_DATA_QUASI(CLEANINPUTS.out, perms.combine(samplenums))
  DESEQ_QUASI(DESEQ_DATA_QUASI.out)
  EDGER_DATA_QUASI(CLEANINPUTS.out, perms.combine(samplenums))
  EDGER_QUASI(EDGER_DATA_QUASI.out)
  // Combine outputs and plot
  QUASI_PLOTS(
    DESEQ_QUASI.out.collect(),
    EDGER_QUASI.out.collect()
  )

}
