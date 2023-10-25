process EDGER_DATA_QUASI{

  input:
  path edger_in
  tuple path(samplesheet), path(countsframe)
  val x

  output:
  tuple path("countsframe_quasi.csv"), path("trueDEGs.RData"), path(samplesheet)

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse", quietly = TRUE)
  library("DESeq2", quietly = TRUE)

  samplesheet = read.csv("$samplesheet")
  countsframe.clean = read.csv("$countsframe", row.names = 1)
  edger.res = read.csv("$edger_in", row.names = 1)
  set.seed($x)

  # Identify 'true' degs with very strict FDR
  truedegs = row.names(subset(edger.res,padj<0.000001))
  # Select half of true degs
  keepnum = round(length(truedegs)/2)

  # Temporary lines for testing 
  DEBUG=FALSE
  if(DEBUG){keepnum = 100 ; truedegs = sample(row.names(results(dds.deg)),200) } 

  #Select a subset of the 'true' degs to keep constant
  keepdegs = sample(truedegs, size = keepnum, replace = FALSE)

  # Permute everything...
  permuterow = function(x){ sample(x, size = ncol(countsframe.clean)) }
  countsframe.quasi = apply(X = countsframe.clean, MARGIN = 1, FUN = permuterow) %>% 
    t() %>% 
    data.frame() %>% 
    `colnames<-`(colnames(countsframe.clean))
  # Then restore true counts for the 'true' DEGs
  countsframe.quasi[keepdegs,] = countsframe.clean[keepdegs,]

  # Output new quasi-permuted frame *and* the set of true degs associated with that frame
  write.csv(countsframe.quasi,"countsframe_quasi.csv")
  save(keepdegs, file="trueDEGs.RData")
  """
}