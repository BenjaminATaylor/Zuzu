process DESEQ_PERMUTE{
  
  label 'mid_job'
  tag "replicate=${x}"

  input: 
  tuple val(x), path(samplesheet), path(countsframe)

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
  nDEGs = nrow(subset(deseq.results,padj<0.05))
  save(nDEGs,file = "nDEGs.RData")
  write.csv(data.frame(deseq.results),row.names=TRUE, file = "deseq_table.csv")
  """
}
