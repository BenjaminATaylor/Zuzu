process PERMUTE_PLOTS{

  debug false

  publishDir "$params.outdir", pattern: 'permute_plot.pdf'
  publishDir "$params.outdir", pattern: 'permute_results.csv', mode: 'copy'

  input:
  path deseq_table
  path(deseq_perms, stageAs: 'deseq_permute/*')
  path edger_table
  path(edger_perms, stageAs: 'edger_permute/*')
  path wilcox_table
  path(wilcox_perms, stageAs: 'wilcox_permute/*')
  path SVC_table, stageAs: 'svc_observed.csv'
  path(SVC_perms, stageAs: 'svc_permute/*')
  path enet_summary, stageAs: 'enet_observed.csv'
  path(enet_perms, stageAs: 'enet_permute/*')

  output:
  path "permute_plot.pdf"
  path "permute_results.csv"

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse", quietly = TRUE)

  # Read one conventional method's keyed permutation tables
  read.permuted.tables = function(directory, method){
    files = list.files(directory, full.names = TRUE)
    if(length(files) == 0L){
      stop(sprintf("No %s permutation tables were supplied", method))
    }

    lapply(files, function(file){
      result = read.csv(file, row.names = 1, check.names = FALSE)
      replicate.match = str_match(
        basename(file),
        "replicate_([0-9]+)_table\\\\.csv\$"
      )[,2]
      if(is.na(replicate.match) || !"padj" %in% names(result)){
        stop(sprintf("Invalid %s permutation table: %s", method, basename(file)))
      }
      data.frame(
        method = method,
        label_assignment = "permuted",
        replicate = as.integer(replicate.match),
        n_tested = nrow(result),
        n_selected = sum(result\$padj < 0.05, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }) %>% bind_rows()
  }

  # Record the observed and permuted results for the conventional methods
  deseq.table = read.csv("$deseq_table", row.names = 1, check.names = FALSE)
  edger.table = read.csv("$edger_table", row.names = 1, check.names = FALSE)
  wilcox.table = read.csv("$wilcox_table", row.names = 1, check.names = FALSE)
  observed.results = bind_rows(
    data.frame(
      method = "DESeq2", label_assignment = "observed", replicate = NA_integer_,
      n_tested = nrow(deseq.table),
      n_selected = sum(deseq.table\$padj < 0.05, na.rm = TRUE)
    ),
    data.frame(
      method = "edgeR", label_assignment = "observed", replicate = NA_integer_,
      n_tested = nrow(edger.table),
      n_selected = sum(edger.table\$padj < 0.05, na.rm = TRUE)
    ),
    data.frame(
      method = "Wilcoxon", label_assignment = "observed", replicate = NA_integer_,
      n_tested = nrow(wilcox.table),
      n_selected = sum(wilcox.table\$padj < 0.05, na.rm = TRUE)
    )
  )
  permuted.results = bind_rows(
    read.permuted.tables("deseq_permute", "DESeq2"),
    read.permuted.tables("edger_permute", "edgeR"),
    read.permuted.tables("wilcox_permute", "Wilcoxon")
  )

  # Include SVC candidates only when specified by the user
  if("$params.mlstep" == "true"){
    SVC.table = read.csv("$SVC_table", row.names = 1, check.names = FALSE)
    observed.results = bind_rows(
      observed.results,
      data.frame(
        method = "SVC RFE-CV candidates",
        label_assignment = "observed",
        replicate = NA_integer_,
        n_tested = nrow(SVC.table),
        n_selected = sum(SVC.table\$selected_feature == 1, na.rm = TRUE)
      )
    )

    svc.files = list.files("svc_permute", full.names = TRUE)
    if(length(svc.files) == 0L){
      stop("No SVC permutation feature counts were supplied")
    }
    SVC.permuted.results = lapply(svc.files, function(file){
      replicate.match = str_match(
        basename(file),
        "replicate_([0-9]+)_n_selected_features\\\\.txt\$"
      )[,2]
      if(is.na(replicate.match)){
        stop(sprintf("Invalid SVC permutation count: %s", basename(file)))
      }
      data.frame(
        method = "SVC RFE-CV candidates",
        label_assignment = "permuted",
        replicate = as.integer(replicate.match),
        n_tested = nrow(SVC.table),
        n_selected = as.integer(read.delim(file, header = FALSE)[1,1])
      )
    }) %>% bind_rows()
    permuted.results = bind_rows(permuted.results, SVC.permuted.results)
  }

  # Include the fixed provisional elastic-net selector when requested
  if("$params.enetstep" == "true"){
    enet.observed = read.csv("$enet_summary", check.names = FALSE)
    stopifnot(
      nrow(enet.observed) == 1L,
      all(c("n_tested", "stable_feature_count") %in% names(enet.observed))
    )
    observed.results = bind_rows(
      observed.results,
      data.frame(
        method = "Stability-selected elastic net",
        label_assignment = "observed",
        replicate = NA_integer_,
        n_tested = as.integer(enet.observed\$n_tested[[1]]),
        n_selected = as.integer(enet.observed\$stable_feature_count[[1]])
      )
    )

    enet.files = list.files("enet_permute", full.names = TRUE)
    if(length(enet.files) == 0L){
      stop("No elastic-net permutation summaries were supplied")
    }
    enet.permuted.results = lapply(enet.files, function(file){
      frame = read.csv(file, check.names = FALSE)
      replicate.match = str_match(
        basename(file),
        "replicate_([0-9]+)_summary\\\\.csv\$"
      )[,2]
      if(is.na(replicate.match) || nrow(frame) != 1L ||
         !all(c("n_tested", "stable_feature_count") %in% names(frame))){
        stop(sprintf("Invalid elastic-net permutation summary: %s", basename(file)))
      }
      data.frame(
        method = "Stability-selected elastic net",
        label_assignment = "permuted",
        replicate = as.integer(replicate.match),
        n_tested = as.integer(frame\$n_tested[[1]]),
        n_selected = as.integer(frame\$stable_feature_count[[1]])
      )
    }) %>% bind_rows()
    permuted.results = bind_rows(permuted.results, enet.permuted.results)
  }

  # Check and retain the complete plot data for later cross-run comparisons
  permute.results = bind_rows(observed.results, permuted.results) %>%
    mutate(selection_proportion = n_selected / n_tested)
  if(any(!is.finite(permute.results\$selection_proportion)) ||
     any(permute.results\$n_selected < 0) ||
     any(permute.results\$n_selected > permute.results\$n_tested)){
    stop("Permutation result counts are inconsistent with their tested universes")
  }
  universe.counts = permute.results %>%
    group_by(method) %>%
    summarize(n_universes = n_distinct(n_tested), .groups = "drop")
  if(any(universe.counts\$n_universes != 1L) ||
     n_distinct(permute.results\$n_tested) != 1L){
    stop("Observed and permuted methods do not share one tested gene universe")
  }
  if(anyDuplicated(
    subset(permute.results, label_assignment == "permuted") %>%
      select(method, replicate)
  )){
    stop("Permutation result table contains duplicate method/replicate keys")
  }

  method.levels = unique(permute.results\$method)
  permute.results = permute.results %>%
    mutate(method = factor(method, levels = method.levels),
           label_assignment = factor(
             label_assignment,
             levels = c("observed", "permuted")
           )) %>%
    arrange(method, label_assignment, replicate)
  write.csv(permute.results, "permute_results.csv", row.names = FALSE)

  # Plot each observed selection count over its permutation distribution
  print(gg.permute <- ggplot(
    subset(permute.results, label_assignment == "permuted"),
    aes(x = method, y = n_selected)
  ) +
    geom_point(size = 3, alpha = 0.7) +
    geom_point(
      data = subset(permute.results, label_assignment == "observed"),
      size = 4, color = "black", fill = "red", shape = 23
    ) +
    labs(x = "Method", y = "Number of selected genes") +
    theme_bw()
  )

  ggsave(gg.permute,
       filename = "permute_plot.pdf",
       device = "pdf", bg = "transparent",
       width =  30, height = 20, units = "cm")

  """

}
