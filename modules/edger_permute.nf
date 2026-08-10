process EDGER_PERMUTE{

  tag "replicate=${x}"

  input: 
  tuple val(x), path(samplesheet), path(countsframe)

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

  # Restore sample sheet factor levels
  reflevel = "$params.reflevel"
  samplesheet\$phenotype = relevel(factor(samplesheet\$phenotype), ref = reflevel)
  stopifnot(nlevels(samplesheet\$phenotype) == 2L,
            identical(levels(samplesheet\$phenotype)[1], reflevel))
  altlevel = setdiff(levels(samplesheet\$phenotype), reflevel)
  stopifnot(length(altlevel) == 1L)

  # Basic edgeR workflow, taken from the documentation
  dds.edge = DGEList(counts=countsframe.clean,group=samplesheet\$phenotype)
  dds.edge = normLibSizes(dds.edge)
  design.edge = model.matrix(~phenotype, data = samplesheet)
  phenotype.coef = which(attr(design.edge, "assign") == 1L)
  stopifnot(length(phenotype.coef) == 1L)
  stopifnot(all(design.edge[samplesheet\$phenotype == reflevel, phenotype.coef] == 0),
            all(design.edge[samplesheet\$phenotype == altlevel, phenotype.coef] == 1))
  edge.contrast = rep(0, ncol(design.edge))
  edge.contrast[phenotype.coef] = 1
  dds.edge.est = estimateDisp(dds.edge,design.edge)
  edge.fit = glmQLFit(dds.edge.est,design.edge)
  edge.qlf = glmQLFTest(edge.fit,contrast=edge.contrast)

  # Tabular output (perform FDR correction directly rather than relying on topTags to do it for us)
  edge.out = edge.qlf\$table %>% mutate(padj = p.adjust(PValue, method = "BH")) 

  nDEGs = nrow(subset(edge.out,padj<0.05))
  save(nDEGs,file = "nDEGs.RData")
  write.csv(edge.out, file="edger_table.csv", row.names = T)
  """
}
