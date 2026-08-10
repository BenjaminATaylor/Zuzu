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

  # Restore sample sheet factor levels
  reflevel = "$params.reflevel"
  samplesheet\$phenotype = relevel(factor(samplesheet\$phenotype), ref = reflevel)
  stopifnot(nlevels(samplesheet\$phenotype) == 2L,
            identical(levels(samplesheet\$phenotype)[1], reflevel))
  altlevel = setdiff(levels(samplesheet\$phenotype), reflevel)
  stopifnot(length(altlevel) == 1L)

  #generate model
  dds = DESeqDataSetFromMatrix(countData = countsframe.clean,
                               colData = samplesheet,
                               design = as.formula(~phenotype))
  # Run the default analysis for DESeq2
  dds.deg = DESeq(dds, fitType = "parametric", betaPrior = FALSE)
  deseq.results = results(dds.deg,
                          contrast = c("phenotype", altlevel, reflevel))
  save(dds.deg, file = "dds_deg.RData")
  write.csv(data.frame(deseq.results),row.names=TRUE, file = "deseq_table.csv")
  """
}
