process DESEQ_PERMUTE{

  input: 
  tuple path(samplesheet), path(countsframe) 

  output:
  path "deseq_table.csv", emit: outfile
  path "nDEGs.RData", emit: nDEGs

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse")
  library("DESeq2")

  samplesheet = read.csv("$samplesheet")
  countsframe.clean = read.csv("$countsframe", row.names = 1, check.names = FALSE)

  #generate model
  dds = DESeqDataSetFromMatrix(countData = countsframe.clean,
                               colData = samplesheet,
                               design = as.formula(~phenotype))
  # Run the default analysis for DESeq2
  dds.deg = DESeq(dds, fitType = "parametric", betaPrior = FALSE)
  nDEGs = nrow(subset(results(dds.deg),padj<0.05))
  save(nDEGs,file = "nDEGs.RData")
  write.csv(data.frame(results(dds.deg)),row.names=TRUE, file = "deseq_table.csv")
  """
}
