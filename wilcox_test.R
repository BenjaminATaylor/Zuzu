setwd("/Users/benjamin/Repositories/Zuzu/work/c2/b60f7cd584f4e4c9805ce3ca83f0c1/")

#!/usr/bin/env Rscript
library("tidyverse", quietly = TRUE)
library("edgeR", quietly = TRUE)

samplesheet = read.csv("samplesheet.csv")
countsframe.clean = read.csv("countsframe_clean.csv", row.names = 1, check.names = F)
set.seed(3)

##Use edgeR by TMM normalization and transfer to CPM (Counts Per Million)
countsframe.dge = DGEList(counts=countsframe.clean,group=samplesheet$phenotype)
countsframe.dge = calcNormFactors(countsframe.dge,method="TMM")
countsframe.dge = cpm(countsframe.dge) %>% as.data.frame()
allgenes = row.names(countsframe.dge)

# Split data into ref and alt counts
refcounts = select(countsframe.dge, subset(samplesheet, phenotype == reflevel)$sample)
altcounts = select(countsframe.dge, subset(samplesheet, phenotype != reflevel)$sample)
# Run wilcox tests on each gene
gene.wilcox = function(x){
  this.wilcox = wilcox.test(as.numeric(refcounts[x,]), as.numeric(altcounts[x,]))
  return(this.wilcox$p.value)
}
wilcox.p = sapply(allgenes, gene.wilcox)
# Adjust to FDR
wilcox.padj = p.adjust(wilcox.p,method = "BH") %>% 
  `names<-`(allgenes) %>% 
  data.frame() %>%
  `colnames<-`("FDR")

write.csv(wilcox.padj, file = "wilcox_out.csv")

