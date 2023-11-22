process PERMUTE_HISTS{

  publishDir "$params.outdir"

  input: 
  val deseq_infiles
  val edger_infiles
  val wilcox_infiles

  output:
  path 'permute_hist.pdf'

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse", quietly = TRUE)
  library("reshape2", quietly = TRUE)

  #print("$deseq_infiles")
  #print("$edger_infiles")
  nperms = as.numeric("$params.nperms")

  readfun = function(x){
    res <- read.csv(x,row.names = 1)
    df <- data.frame(result = (res\$padj<0.05), row.names = row.names(res))
    return(df)
  }

  deseq.inlist = str_remove_all("$deseq_infiles","[\\\\[\\\\] ]") %>% 
    strsplit(split = ",") %>% unlist()
  deseq.permdegs = lapply(deseq.inlist,readfun) %>% 
    list_cbind() %>%
    rowSums()

  edger.inlist = str_remove_all("$edger_infiles","[\\\\[\\\\] ]") %>% 
    strsplit(split = ",") %>% unlist()
  edger.permdegs = lapply(edger.inlist,readfun) %>% 
    list_cbind() %>%
    rowSums()

  wilcox.inlist = str_remove_all("$wilcox_infiles","[\\\\[\\\\] ]") %>% 
    strsplit(split = ",") %>% unlist()
  wilcox.permdegs = lapply(wilcox.inlist,readfun) %>% 
    list_cbind() %>%
    rowSums()

  permhist.in = cbind(deseq.permdegs,edger.permdegs,wilcox.permdegs) %>%
    `colnames<-`(c("DESeq2","edgeR","Wilcoxon")) %>%
    melt()

  # Generate breaks contingent on number of reps
  cutbreaks = cut(seq(1,nperms+1),breaks = 5, right = TRUE) %>% levels() %>%
    str_match("\\\\([^,]*") %>%
    str_remove("\\\\(") %>%
    as.numeric() %>%
    round()
  #deal with an idiosyncrasy that causes cut to return 0 as the start value of the first bin instead of 1
  cutbreaks[1] = 1
  
  #plot
  gg.permhist = ggplot(permhist.in, aes(x = value)) +
    geom_histogram(color = "black", fill = "white", bins = 100) +
    scale_x_continuous(limits = c(1,(nperms+1)), breaks = cutbreaks) +
    labs(x = "Number of permutations in which gene was\nincorrectly identified as a DEG",
         y = "Number of genes") +
    theme_bw() +
    facet_grid(~Var2) +
    theme(strip.background = element_blank())

  ggsave(gg.permhist, 
   filename = "permute_hist.pdf",
   device = "pdf", bg = "transparent",
   width =  30, height = 20, units = "cm")

  """
}