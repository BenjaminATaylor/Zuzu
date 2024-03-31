setwd("/Users/benjamin/Repositories/Zuzu/work/0c/e9177ffc36dc9da23286d9eae6d01e")

#!/usr/bin/env Rscript
library("tidyverse", quietly = TRUE)
library("RedRibbon", quietly = TRUE)
library("UpSetR", quietly = TRUE)

#hypergeometric overlap test function
gene_overlap_test = function(list1, list2, background,verbose=TRUE,lower.tail=FALSE){
  n_A = length(list1)
  n_B = length(list2)
  n_C = length(background)
  n_A_B = length(intersect(list1,list2))
  hyp = phyper(n_A_B - 1, n_A, n_C-n_A, n_B, lower.tail = lower.tail)
  expc = round((n_A/n_C)*n_B)
  jac = n_A_B/(n_A+n_B-n_A_B)
  if(verbose){print(paste0(n_A_B," matching from lists of lengths ",n_A," and ",
                           n_B,"; p=",round(hyp,4),". Expected intersection length = ",expc))}
  res = list(hypergeom = round(hyp,4),
             jaccard = jac,
             intersect = n_A_B)
  return(res)
}

# Get DEGs from full model
deseq.table = read.csv("deseq_table.csv", row.names = 1)
edger.table = read.csv("edger_table.csv", row.names = 1)
wilcox.table = read.csv("wilcox_table.csv", row.names = 1)
getdegs = function(df, ml=FALSE){ 
  if(ml){
    return(row.names(subset(df, DEGstatus==1))) 
  } else { 
    return(row.names(subset(df, padj<0.05))) 
  } 
}
degs.deseq = getdegs(deseq.table)
degs.edger = getdegs(edger.table)
degs.wilcox = getdegs(wilcox.table)
# Only include ML outputs if specified by user
# TODO: Fix this to work for ML outputs also
if("false" == "true"){
  svc.table = read.csv("svc_table.csv", row.names = 1)
  degs.svc = getdegs("svc_table.csv", ml=TRUE)
}

compframe = cbind(select(deseq.table,padj), 
                  select(edger.table,padj),
                  select(wilcox.table,padj))  %>% 
  remove_rownames() %>%
  `colnames<-`(c("deseq","edger","wilcox"))

pairs = t(combn(colnames(compframe),2))

for(i in 1:nrow(pairs)){
  
  a = compframe[,pairs[i,1]]
  b = compframe[,pairs[i,2]]
  
  # First let's test the correlations of the p-values of the two sets (only relevant where we have pvals)
  # RRHO would be nice here but neither of the two existing packages (RRHO and RedRibbon) are suitable- the former due to bugs and the latter due to instability
  # TODO: apply my own RRHO algorithm
  cor.kendall = cor.test(a, b, 
                         alternative = "two.sided",
                         method = "kendall")$estimate
  cor.pearson = cor.test(a, b, 
                         alternative = "two.sided",
                         method = "pearson")$estimate
  plot(a,b)
  
  # As expected, the pval correlations themselves are usually very high. The more interesting question is how similar each pair of methods is in terms of the specific set of DEGs identified
 
  alldegs = row.names(compframe)
  a = get(paste0("degs.",pairs[i,1]))
  b = get(paste0("degs.",pairs[i,2]))
  
  thisgot = gene_overlap_test(a, b, alldegs, verbose = FALSE)
  jacc = thisgot$jaccard
  hyprp = thisgot$hypergeom
  
}

# also generate an upSet showing the intersections
listInput.combined <- list("DESeq2" = degs.deseq,
                           "edgeR" = degs.edger,
                           "Wilcoxon" = degs.wilcox)

print(upset.combined <- upset(fromList(listInput.combined), order.by = "degree", 
                              point.size = 4, line.size = 2, 
                              mainbar.y.label = "DEG intersections", 
                              sets.x.label = "Number of DEGs", 
                              empty.intersections = "on", text.scale = 3.2)
)

pdf(file="upSet_basic.pdf",width = 24, height = 18) ; upset.combined ;dev.off()







