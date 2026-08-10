process EDGER_DATA_QUASI{

  input:
  tuple path(samplesheet), path(countsframe)
  tuple val(x), val(refnum), val(altnum)

  output:
  tuple path("countsframe_quasi.csv"),
        path("trueDEGs.RData"),
        path("samplesheet_quasi.csv"),
        val(x),
        val(refnum),
        val(altnum)

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse", quietly = TRUE)
  library("edgeR", quietly = TRUE)

  samplesheet = read.csv("$samplesheet")
  countsframe.clean = read.csv("$countsframe", row.names = 1, check.names = F)
  set.seed($x)
  ref.num = as.integer($refnum)
  alt.num = as.integer($altnum)

  # Take a random subset of the samples based on our chosen subsample size
  phenos = relevel(factor(unique(samplesheet\$phenotype)), ref = "$params.reflevel")
  chosen.samples = c(sample(subset(samplesheet, phenotype == levels(phenos)[1])\$sample,ref.num),
                    sample(subset(samplesheet, phenotype == levels(phenos)[2])\$sample,alt.num))
  samplesheet.sub = subset(samplesheet, sample %in% chosen.samples)
  countsframe.sub = countsframe.clean[,samplesheet.sub\$sample]
  stopifnot(nrow(samplesheet.sub) == ref.num + alt.num)
  stopifnot(sum(samplesheet.sub\$phenotype == "$params.reflevel") == ref.num)
  stopifnot(sum(samplesheet.sub\$phenotype != "$params.reflevel") == alt.num)
  stopifnot(identical(colnames(countsframe.sub), samplesheet.sub\$sample))

  # Re-run edgeR for our new dataset
  dds.edge = DGEList(counts=countsframe.sub,group=samplesheet.sub\$phenotype)
  dds.edge = normLibSizes(dds.edge)
  design.edge = model.matrix(~samplesheet.sub\$phenotype)
  dds.edge.est = estimateDisp(dds.edge,design.edge)
  edge.fit = glmQLFit(dds.edge.est,design.edge)
  edge.qlf = glmQLFTest(edge.fit,coef=2)
  edger.res = edge.qlf\$table %>% mutate(padj = p.adjust(PValue, method = "BH")) 

  # Identify 'true' degs with very strict FDR
  truedegs = row.names(subset(edger.res,padj<$params.truecutoff))
  # Select half of true degs
  keepnum = round(length(truedegs)/2)

  # Temporary lines for testing 
  DEBUG=FALSE
  if(DEBUG){keepnum = 100 ; truedegs = sample(row.names(results(dds.deg)),200) } 

  #Select a subset of the 'true' degs to keep constant
  keepdegs = sample(truedegs, size = keepnum, replace = FALSE)

  # Permute only the selected samples so anchor discovery and recovery use the
  # same dataset size.
  permuterow = function(x){ sample(x, size = ncol(countsframe.sub)) }
  countsframe.quasi = apply(X = countsframe.sub, MARGIN = 1, FUN = permuterow) %>%
    t() %>% 
    data.frame() %>% 
    `colnames<-`(colnames(countsframe.sub))
  # Then restore the selected-sample counts for the retained anchors.
  countsframe.quasi[keepdegs,] = countsframe.sub[keepdegs,]
  stopifnot(identical(colnames(countsframe.quasi), samplesheet.sub\$sample))

  # Output the quasi-permuted frame, its anchors, and the matching sample subset.
  write.csv(countsframe.quasi,"countsframe_quasi.csv")
  write.csv(samplesheet.sub,"samplesheet_quasi.csv", row.names = FALSE)
  save(keepdegs, file="trueDEGs.RData")
  """
}
