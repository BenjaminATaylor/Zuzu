process QUASI_PLOTS{

  publishDir "$params.outdir"

  input:
  path(deseq_quasis,  stageAs: 'deseq/outframe_?.csv')
  path(edger_quasis,  stageAs: 'edger/outframe_?.csv')
  path(wilcox_quasis, stageAs: 'wilcox/outframe_?.csv')
  path(svm_quasis,    stageAs: 'svc/outframe_?.csv')

  output:
  path "quasi_plot.pdf"
  path "quasi_results.csv"

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse", quietly = TRUE)
  library("reshape2", quietly = TRUE)

  read_quasis = function(dir, label){
    lapply(list.files(dir, full.names = TRUE), read.csv) %>%
      bind_rows() %>%
      mutate(method = label)
  }

  # Combine method-specific replicates into one ordered results table
  quasi.results = bind_rows(
    read_quasis("deseq",  "DESeq2"),
    read_quasis("edger",  "edgeR"),
    read_quasis("wilcox", "Wilcoxon")
  ) %>%
    arrange(method, replicate, totalnum, refnum, altnum)

  write.csv(quasi.results, "quasi_results.csv", row.names = FALSE)

  # Reshape the self-consistency metrics for faceted plotting
  gg.quasi.input = melt(
    quasi.results,
    id.vars = c("method", "replicate", "refnum", "altnum", "totalnum"),
    measure.vars = c("power", "FDP", "youden_index")
  ) %>%
    mutate(variable = fct_recode(variable,
                                 "Youden's index" = "youden_index")) %>%
    mutate(sample_size = paste0(refnum, " ref + ", altnum, " alt")) %>%
    arrange(totalnum, refnum, altnum) %>%
    mutate(sample_size = factor(sample_size, levels = unique(sample_size)))

  gg.quasi = ggplot(gg.quasi.input, aes(x = method, y = value)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_y_continuous(limits = c(-1,1)) +
    labs(x = "Method", y = "Self-consistency metric") +
    facet_grid(sample_size~variable) +
    theme_bw() +
    theme(strip.background = element_blank())

  ggsave(gg.quasi,
     filename = "quasi_plot.pdf",
     device = "pdf", bg = "transparent",
     width = 35, height = 20, units = "cm")
  """
}
