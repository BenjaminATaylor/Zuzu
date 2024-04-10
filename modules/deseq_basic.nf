process DESEQ_BASIC{
  
  label 'mid_job'

  input: 
  tuple path(samplesheet), path(countsframe) 

  output:
  path "dds_deg.RData", emit: dds
  path "deseq_table.csv", emit: table

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse", quietly = TRUE)
  library("DESeq2", quietly = TRUE)

  samplesheet = read.csv("$samplesheet")
  countsframe.clean = read.csv("$countsframe", row.names = 1)

  #generate model
  dds = DESeqDataSetFromMatrix(countData = countsframe.clean,
                               colData = samplesheet,
                               design = as.formula(~phenotype))
  # Run the default analysis for DESeq2
  dds.deg = DESeq(dds, fitType = "parametric", betaPrior = FALSE)
  save(dds.deg, file = "dds_deg.RData")
  write.csv(data.frame(results(dds.deg)),row.names=TRUE, file = "deseq_table.csv")
  """
}