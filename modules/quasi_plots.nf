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
    id.vars = c("method", "replicate", "refnum", "altnum", "totalnum", "n_truth"),
    measure.vars = c("power", "FDP", "youden_index")
  ) %>%
    mutate(variable = fct_recode(variable,
                                 "Youden's index" = "youden_index")) %>%
    mutate(sample_size = paste0(refnum, " ref + ", altnum, " alt")) %>%
    arrange(totalnum, refnum, altnum) %>%
    mutate(sample_size = factor(sample_size, levels = unique(sample_size))) %>%
    mutate(anchor_count = case_when(
      n_truth == 0  ~ "No anchors",
      n_truth <= 10 ~ "1-10 anchors",
      TRUE          ~ ">10 anchors"
    )) %>%
    mutate(anchor_count = factor(
      anchor_count,
      levels = c(">10 anchors", "1-10 anchors", "No anchors")
    ))

  # Mark results calculated from small or empty anchor sets
  gg.quasi = ggplot(gg.quasi.input, aes(x = method, y = value, shape = anchor_count)) +
    geom_point(size = 3, alpha = 0.7, na.rm = TRUE) +
    scale_shape_manual(values = c(
      ">10 anchors" = 16,
      "1-10 anchors" = 1,
      "No anchors" = 4
    )) +
    scale_y_continuous(limits = c(-1,1)) +
    labs(x = "Method",
         y = "Self-consistency metric",
         shape = "Retained anchors",
         caption = paste(
           "Power and Youden's index are undefined when no anchors are retained;",
           "Youden's index is also undefined when no background genes remain."
         )) +
    facet_grid(sample_size~variable) +
    theme_bw() +
    theme(strip.background = element_blank())

  ggsave(gg.quasi,
     filename = "quasi_plot.pdf",
     device = "pdf", bg = "transparent",
     width = 35, height = 20, units = "cm")
  """
}
