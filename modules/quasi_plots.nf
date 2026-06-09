process QUASI_PLOTS{

  publishDir "$params.outdir"

  input:
  path(deseq_quasis,  stageAs: 'deseq/outframe_?.csv')
  path(edger_quasis,  stageAs: 'edger/outframe_?.csv')
  path(wilcox_quasis, stageAs: 'wilcox/outframe_?.csv')
  path(svm_quasis,    stageAs: 'svc/outframe_?.csv')

  output:
  path "quasi_plot.pdf"

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse", quietly = TRUE)
  library("reshape2", quietly = TRUE)

  read_quasis = function(dir, label){
    lapply(list.files(dir, full.names = TRUE), read.csv) %>%
      bind_rows() %>%
      mutate_all(as.numeric) %>%
      mutate(method = label) %>%
      melt(id.vars = c("method", "samplenum"))
  }

  gg.quasi.input = rbind(
    read_quasis("deseq",  "DESeq2"),
    read_quasis("edger",  "edgeR"),
    read_quasis("wilcox", "Wilcoxon")
  )

  gg.quasi = ggplot(gg.quasi.input, aes(x = method, y = value)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_y_continuous(limits = c(0,1)) +
    labs(x = "Method", y = "Replicates per group") +
    facet_grid(samplenum~variable) +
    theme_bw() +
    theme(strip.background = element_blank())

  ggsave(gg.quasi,
     filename = "quasi_plot.pdf",
     device = "pdf", bg = "transparent",
     width = 35, height = 20, units = "cm")
  """
}
