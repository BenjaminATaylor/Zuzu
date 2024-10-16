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