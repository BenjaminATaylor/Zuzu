process DESEQ_DATA_QUASI{

  input:
  tuple path(samplesheet), path(countsframe)
  tuple val(x), val(samplenum)

  output:
  tuple path("countsframe_quasi.csv"), path("trueDEGs.RData"), path(samplesheet), val(samplenum)

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse", quietly = TRUE)
  library("DESeq2", quietly = TRUE)

  samplesheet = read.csv("$samplesheet")
  countsframe.clean = read.csv("$countsframe", row.names = 1, check.names = F)
  set.seed($x)
  n.samples = $samplenum

  # Take a random subset of the samples based on our chosen subsample size
  phenos = unique(samplesheet\$phenotype)
  chosen.samples = c(sample(subset(samplesheet, phenotype == phenos[1])\$sample,n.samples),
                     sample(subset(samplesheet, phenotype == phenos[2])\$sample,n.samples))
  samplesheet.sub = subset(samplesheet, sample %in% chosen.samples)
  countsframe.sub = countsframe.clean[,samplesheet.sub\$sample]
  stopifnot(colnames(countsframe.sub) == samplesheet.sub\$sample)

  # Re-run deseq for our new dataset
  dds = DESeqDataSetFromMatrix(countData = countsframe.sub,
                               colData = samplesheet.sub,
                               design = as.formula(~phenotype))
  # Run the default analysis for DESeq2
  dds.deg = DESeq(dds, fitType = "parametric", betaPrior = FALSE)

  # Identify 'true' degs with very strict FDR
  truedegs = row.names(subset(results(dds.deg,alpha = 0.000001),padj<0.000001))
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