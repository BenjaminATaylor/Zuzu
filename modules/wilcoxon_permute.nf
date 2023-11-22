process WILCOXON_PERMUTE{

  input: 
  tuple path(samplesheet), path(countsframe) 

  output:
  path "wilcox_table.csv", emit: outfile
  path "nDEGs.RData", emit: nDEGs

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse")
  library("edgeR")

  samplesheet = read.csv("$samplesheet")
  countsframe.clean = read.csv("$countsframe", row.names = 1, check.names = FALSE)
  reflevel = "$params.reflevel"

  ##Use edgeR by TMM normalization and transfer to CPM (Counts Per Million)
  countsframe.dge = DGEList(counts=countsframe.clean,group=samplesheet\$phenotype)
  countsframe.dge = calcNormFactors(countsframe.dge,method="TMM")
  countsframe.dge = cpm(countsframe.dge) %>% as.data.frame()
  allgenes = row.names(countsframe.dge)

  # Split data into ref and alt counts
  refcounts = select(countsframe.dge, subset(samplesheet, phenotype == reflevel)\$sample)
  altcounts = select(countsframe.dge, subset(samplesheet, phenotype != reflevel)\$sample)
  # Run wilcox tests on each gene
  gene.wilcox = function(x){
    this.wilcox = wilcox.test(as.numeric(refcounts[x,]), as.numeric(altcounts[x,]))
    return(this.wilcox\$p.value)
  }
  wilcox.p = sapply(allgenes, gene.wilcox)
  # Adjust to FDR
  wilcox.padj = p.adjust(wilcox.p,method = "BH") %>% 
    `names<-`(allgenes) %>% 
    data.frame() %>%
    `colnames<-`("padj")

  nDEGs = length(which(wilcox.padj<0.05))
  save(nDEGs,file = "nDEGs.RData")
  write.csv(wilcox.padj, row.names = TRUE, file = "wilcox_table.csv")
  """
}

