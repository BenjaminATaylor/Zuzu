#!/usr/bin/env nextflow

// Modules to include
include { CLEANINPUTS } from './modules/cleaninputs.nf'
include { QUALITYCONTROL } from './modules/qualitycontrol.nf'
include { DESEQ_BASIC } from './modules/deseq_basic.nf'
include { DATA_PERMUTE } from './modules/data_permute.nf'
include { DESEQ_PERMUTE } from './modules/deseq_permute.nf'
include { DATA_QUASI } from './modules/data_quasi.nf'
include { DESEQ_QUASI } from './modules/deseq_quasi.nf'
include { DESEQ_QUASI_COLLECT } from './modules/deseq_quasi_collect.nf'

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
  path "edger_table.csv"
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

process PERMUTE_PLOTS{

  debug true
  publishDir "$params.outdir"

  input:
  path deseq_table
  val deseq_perms
  path edger_table
  val edger_perms

  output:
  path "permute_plot.pdf"

  """
  #!/usr/bin/env Rscript
  library("tidyverse")
  library("DESeq2")

  ## DESeq inputs
  # DEGs from full model
  deseq.res = get(load("$deseq_table"))
  deseq.table = results(deseq.res)
  deseq.ndegs = nrow(subset(deseq.table, padj<0.05))
  # Permutations
  deseq.inlist = str_remove_all("$deseq_perms","[\\\\[\\\\] ]") %>% 
    strsplit(split = ",") %>% unlist()
  deseq.perms = sapply(deseq.inlist, function(x) get(load(x))) %>% unname()
  # Collated input
  deseq.permute.input = data.frame( method = "DESeq2", 
                                  nDEGs = deseq.ndegs, 
                                  permutes = deseq.perms
                                  )

  ## edgeR inputs
  # DEGs from full model
  edger.table = read.csv("$edger_table", row.names = 1)
  edger.ndegs = nrow(subset(edger.table, padj<0.05))
  # Permutations
  edger.inlist = str_remove_all("$edger_perms","[\\\\[\\\\] ]") %>% 
    strsplit(split = ",") %>% unlist()
  edger.perms = sapply(edger.inlist, function(x) get(load(x))) %>% unname()
  # Collated input
  edger.permute.input =  data.frame(method = "edgeR", 
                                    nDEGs = edger.ndegs, 
                                    permutes = edger.perms)

  # Combined inputs for plotting
  permuteplot.input = rbind(deseq.permute.input,edger.permute.input)

  print(gg.permute <- ggplot(permuteplot.input, aes(x = method, y = permutes)) +
          geom_point(size = 3, alpha = 0.7) +
          geom_point(aes(x = method, y = nDEGs), 
                    size =4 , color = "black", fill = "red", shape = 23) +
          labs(x = "Method", y = "Number of DEGs") +
          theme_bw()
  )

  ggsave(gg.permute, 
       filename = "permute_plot.pdf",
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
    DESEQ_BASIC.out,
    DESEQ_PERMUTE.out.nDEGs.collect(),
    EDGER_BASIC.out,
    EDGER_PERMUTE.out.nDEGs.collect()
  )

  DATA_QUASI(DESEQ_BASIC.out, CLEANINPUTS.out, perms) |
  DESEQ_QUASI 

  DESEQ_QUASI.out
  .collect() |
  DESEQ_QUASI_COLLECT

  //DESEQ_QUASI_COLLECT.out.view()

}
