process QUASI_PLOTS{

  publishDir "$params.outdir"

  input:
  val deseq_quasis
  val edger_quasis
  val wilcox_quasis
  val svm_quasis

  output:
  path "quasi_plot.pdf"

  """
  #!/usr/bin/env Rscript
  library("tidyverse", quietly = TRUE)
  library("reshape2", quietly = TRUE)

  deseq.inlist = str_remove_all("$deseq_quasis","[\\\\[\\\\] ]") %>% 
    strsplit(split = ",") %>% unlist()

  deseq.quasi.input = 
    sapply(deseq.inlist, read.csv) %>% 
    t() %>% 
    `row.names<-`(NULL) %>% 
    data.frame() %>%
    mutate_all(as.numeric) %>%
    mutate(method = "DESeq2") %>% 
    melt(id.vars = c("method","samplenum"))

  edger.inlist = str_remove_all("$edger_quasis","[\\\\[\\\\] ]") %>% 
    strsplit(split = ",") %>% unlist()

  edger.quasi.input = 
    sapply(edger.inlist, read.csv) %>% 
    t() %>% 
    `row.names<-`(NULL) %>% 
    data.frame() %>%
    mutate_all(as.numeric) %>%
    mutate(method = "edgeR") %>% 
    melt(id.vars = c("method","samplenum"))

  wilcox.inlist = str_remove_all("$wilcox_quasis","[\\\\[\\\\] ]") %>% 
    strsplit(split = ",") %>% unlist()

  wilcox.quasi.input = 
    sapply(wilcox.inlist, read.csv) %>% 
    t() %>% 
    `row.names<-`(NULL) %>% 
    data.frame() %>%
    mutate_all(as.numeric) %>%
    mutate(method = "wilcox") %>% 
    melt(id.vars = c("method","samplenum"))
    
  svm.inlist = str_remove_all("$svm_quasis","[\\\\[\\\\] ]") %>% 
    strsplit(split = ",") %>% unlist()

  svm.quasi.input = 
    sapply(svm.inlist, read.csv) %>% 
    t() %>% 
    `row.names<-`(NULL) %>% 
    data.frame() %>%
    mutate_all(as.numeric) %>%
    mutate(method = "svm") %>% 
    melt(id.vars = c("method","samplenum"))

  gg.quasi.input = rbind(deseq.quasi.input, 
                          edger.quasi.input,
                          wilcox.quasi.input,
                          svm.quasi.input)

  gg.quasi = ggplot(gg.quasi.input, aes(x = method, y = value)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_y_continuous(limits = c(0,1)) +
    labs(x = "Method",y = "Replicates per group") +
    facet_grid(samplenum~variable) +
    theme_bw() +
    theme(strip.background = element_blank())

  ggsave(gg.quasi, 
     filename = "quasi_plot.pdf",
     device = "pdf", bg = "transparent",
     width =  35, height = 20, units = "cm")
  """
}