process FULLSYNTH_PLOTS {

  publishDir "$params.outdir"

  input:
  val deseq_synths
  val edger_synths


  output:
  path "synth_plot.pdf"

  """
  #!/usr/bin/env Rscript
  library("tidyverse", quietly = TRUE)
  library("reshape2", quietly = TRUE)

  deseq.inlist = str_remove_all("$deseq_synths","[\\\\[\\\\] ]") %>% 
    strsplit(split = ",") %>% unlist()

  deseq.synth.input = 
    sapply(deseq.inlist, read.csv) %>% 
    t() %>% 
    `row.names<-`(NULL) %>% 
    data.frame() %>%
    mutate_all(as.numeric) %>%
    mutate(method = "DESeq2") %>% 
    melt(id.vars = c("method","samplenum"))

  edger.inlist = str_remove_all("$edger_synths","[\\\\[\\\\] ]") %>% 
    strsplit(split = ",") %>% unlist()

  edger.synth.input = 
    sapply(edger.inlist, read.csv) %>% 
    t() %>% 
    `row.names<-`(NULL) %>% 
    data.frame() %>%
    mutate_all(as.numeric) %>%
    mutate(method = "edgeR") %>% 
    melt(id.vars = c("method","samplenum"))

    #  wilcox.inlist = str_remove_all("wilcox_synths","[\\\\[\\\\] ]") %>% 
    #    strsplit(split = ",") %>% unlist()
    #
    #  wilcox.synth.input = 
    #    sapply(wilcox.inlist, read.csv) %>% 
    #    t() %>% 
    #    `row.names<-`(NULL) %>% 
    #    data.frame() %>%
    #    mutate_all(as.numeric) %>%
    #    mutate(method = "wilcox") %>% 
    #    melt(id.vars = c("method","samplenum"))

  gg.synth.input = rbind(deseq.synth.input, 
                          edger.synth.input)

  gg.synth = ggplot(gg.synth.input, aes(x = method, y = value)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_y_continuous(limits = c(0,1)) +
    labs(x = "Method",y = "Replicates per group") +
    facet_grid(samplenum~variable) +
    theme_bw() +
    theme(strip.background = element_blank())

  ggsave(gg.synth, 
     filename = "synth_plot.pdf",
     device = "pdf", bg = "transparent",
     width =  30, height = 20, units = "cm")
  """
}
