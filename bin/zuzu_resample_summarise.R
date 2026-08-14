#!/usr/bin/env Rscript

# Parse the deliberately small named-pair command line interface
parse.arguments = function(arguments){
  expected = c("summaries", "genes", "allocations", "libraries", "baseline")
  if(length(arguments) != 2L * length(expected)){
    stop("All resampling summary inputs must be supplied as named pairs")
  }
  keys = sub("^--", "", arguments[seq.int(1L, length(arguments), by = 2L)])
  values = arguments[seq.int(2L, length(arguments), by = 2L)]
  if(!setequal(keys, expected) || anyDuplicated(keys)){
    stop("Unexpected, missing, or duplicated resampling summary argument")
  }
  stats::setNames(as.list(values), keys)
}

arguments = parse.arguments(commandArgs(trailingOnly = TRUE))

# Read every per-dataset table produced by the shared resampling runner
read.directory = function(directory){
  files = list.files(directory, pattern = "[.]csv$", full.names = TRUE)
  if(length(files) == 0L) stop("No CSV files were supplied in ", directory)
  dplyr::bind_rows(lapply(
    files,
    utils::read.csv,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = "character"
  ))
}

suppressPackageStartupMessages(library("tidyverse"))
summaries = read.directory(arguments[["summaries"]])
genes = read.directory(arguments[["genes"]])
allocations = read.directory(arguments[["allocations"]])
libraries = read.directory(arguments[["libraries"]])
baseline = utils::read.csv(
  arguments[["baseline"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)

numeric.summary.columns = c(
  "repeat", "target", "generation_seed", "samples_per_group", "target_depth",
  "reference_n", "alternative_n", "total_n", "n_input_genes", "n_tested"
)
summaries = summaries %>%
  mutate(across(all_of(numeric.summary.columns), as.numeric))
genes = genes %>%
  mutate(across(all_of(c("repeat", "target")), as.numeric))
allocations = allocations %>%
  mutate(across(any_of(c("repeat", "target")), as.numeric))
libraries = libraries %>%
  mutate(across(
    any_of(c(
      "repeat", "target", "original_depth", "target_depth",
      "thinned_prefilter_depth", "postfilter_depth"
    )),
    as.numeric
  ))

# Confirm the compact common schema before combining method results
dataset.keys = c("resampling_type", "repeat", "target")
summary.columns = c(
  dataset.keys, "generation_seed", "samples_per_group", "target_depth",
  "reference_n", "alternative_n", "total_n", "n_input_genes", "n_tested"
)
gene.columns = c(dataset.keys, "gene", "tested")
method.columns = c("DESeq2", "edgeR", "Wilcoxon", "Elastic net")
if(!all(summary.columns %in% names(summaries)) ||
   !all(c(gene.columns, method.columns) %in% names(genes)) ||
   !all(c("gene", method.columns) %in% names(baseline))){
  stop("A resampling input table does not follow the common schema")
}
if(anyDuplicated(summaries[dataset.keys])){
  stop("Resampling dataset keys are not unique")
}

# Convert the complete raw-gene universe to one row per method and dataset
enabled.method.columns = method.columns[vapply(
  genes[method.columns],
  function(values) any(!is.na(values)),
  logical(1)
)]
genes.long = genes %>%
  pivot_longer(
    cols = all_of(enabled.method.columns),
    names_to = "method",
    values_to = "selected",
    values_transform = as.logical
  ) %>%
  filter(!is.na(selected)) %>%
  mutate(tested = as.logical(tested))

baseline.long = baseline %>%
  select(gene, all_of(method.columns)) %>%
  pivot_longer(
    cols = all_of(method.columns),
    names_to = "method",
    values_to = "baseline_selected",
    values_transform = as.logical
  )

genes.long = genes.long %>%
  left_join(
    baseline.long %>% mutate(baseline_tested = TRUE),
    by = c("gene", "method")
  ) %>%
  mutate(
    baseline_tested = replace_na(baseline_tested, FALSE),
    baseline_selected = replace_na(baseline_selected, FALSE)
  )

# Compare every perturbed list with that method's full-data baseline list
results = genes.long %>%
  group_by(across(all_of(c(dataset.keys, "method")))) %>%
  summarise(
    n_selected = sum(selected),
    baseline_n_tested = sum(baseline_tested),
    baseline_n_selected = sum(baseline_selected),
    n_baseline_selected_tested = sum(baseline_selected & tested),
    n_baseline_selected_filtered_out = sum(baseline_selected & !tested),
    intersection = sum(selected & baseline_selected),
    union = sum(selected | baseline_selected),
    n_new_calls = sum(selected & !baseline_selected),
    .groups = "drop"
  ) %>%
  left_join(summaries[,summary.columns], by = dataset.keys) %>%
  mutate(
    selection_proportion = n_selected / n_tested,
    jaccard = if_else(union == 0L, NA_real_, intersection / union),
    baseline_recovery = if_else(
      baseline_n_selected == 0L,
      NA_real_,
      intersection / baseline_n_selected
    ),
    new_call_proportion = if_else(
      n_selected == 0L,
      NA_real_,
      n_new_calls / n_selected
    ),
    empty_list = n_selected == 0L,
    status = case_when(
      n_selected == 0L & baseline_n_selected == 0L ~ "both_lists_empty",
      n_selected == 0L ~ "resample_list_empty",
      baseline_n_selected == 0L ~ "baseline_list_empty",
      TRUE ~ "ok"
    )
  ) %>%
  select(
    resampling_type, `repeat`, generation_seed, target, samples_per_group,
    target_depth, reference_n, alternative_n, total_n, method,
    n_input_genes, n_tested, n_selected, selection_proportion,
    baseline_n_tested, baseline_n_selected, n_baseline_selected_tested,
    n_baseline_selected_filtered_out, intersection, union, jaccard,
    baseline_recovery, n_new_calls, new_call_proportion, empty_list, status
  ) %>%
  mutate(method = factor(method, levels = method.columns)) %>%
  arrange(resampling_type, target, `repeat`, method)

# Retain unconditional and filter-conditional selection frequencies for every gene
gene.frequencies = genes.long %>%
  group_by(resampling_type, target, method, gene) %>%
  summarise(
    baseline_tested = first(baseline_tested),
    baseline_selected = first(baseline_selected),
    n_repeats = n(),
    n_tested_repeats = sum(tested),
    tested_frequency = mean(tested),
    n_selected_repeats = sum(selected),
    selection_frequency = mean(selected),
    selection_frequency_when_tested = if(
      sum(tested) == 0L
    ) NA_real_ else mean(selected[tested]),
    .groups = "drop"
  ) %>%
  mutate(
    samples_per_group = if_else(
      resampling_type == "sample", as.numeric(target), NA_real_
    ),
    target_depth = if_else(
      resampling_type == "depth", as.numeric(target), NA_real_
    )
  ) %>%
  select(
    resampling_type, target, samples_per_group, target_depth, method, gene,
    baseline_tested, baseline_selected, n_repeats, n_tested_repeats,
    tested_frequency, n_selected_repeats, selection_frequency,
    selection_frequency_when_tested
  ) %>%
  mutate(method = factor(method, levels = method.columns)) %>%
  arrange(resampling_type, target, method, gene)

# Write separate copied tables and one compact diagnostic plot for each enabled branch
for(resampling.type in unique(results$resampling_type)){
  type.results = filter(results, resampling_type == resampling.type)
  type.frequencies = filter(
    gene.frequencies,
    resampling_type == resampling.type
  )
  utils::write.csv(
    type.results,
    paste0(resampling.type, "_subsample_results.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    type.frequencies,
    paste0(resampling.type, "_subsample_gene_frequencies.csv"),
    row.names = FALSE
  )

  if(resampling.type == "sample"){
    utils::write.csv(
      filter(allocations, resampling_type == "sample") %>%
        arrange(target, `repeat`, phenotype, sample),
      "sample_subsample_allocations.csv",
      row.names = FALSE
    )
  }
  if(resampling.type == "depth"){
    utils::write.csv(
      filter(libraries, resampling_type == "depth") %>%
        arrange(target, `repeat`, sample),
      "depth_subsample_library_metrics.csv",
      row.names = FALSE
    )
  }

  plot.frame = type.results %>%
    select(method, target, `repeat`, baseline_recovery, jaccard, n_selected) %>%
    pivot_longer(
      cols = c(baseline_recovery, jaccard, n_selected),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(metric = recode(
      metric,
      baseline_recovery = "Baseline-list recovery",
      jaccard = "Jaccard overlap",
      n_selected = "Selected genes"
    ))

  resample.plot = ggplot(
    plot.frame,
    aes(x = target, y = value, group = interaction(method, target))
  ) +
    geom_boxplot(outlier.shape = NA, width = 0.16) +
    geom_point(position = position_jitter(width = 0.04, height = 0), size = 1.2) +
    facet_grid(metric ~ method, scales = "free_y") +
    labs(
      x = if(resampling.type == "sample"){
        "Samples per phenotype"
      } else {
        "Pre-filter count total per sample"
      },
      y = NULL
    ) +
    theme_bw()
  if(resampling.type == "depth"){
    resample.plot = resample.plot + scale_x_log10()
  }
  ggsave(
    paste0(resampling.type, "_subsample_plot.pdf"),
    resample.plot,
    width = 4 * length(unique(plot.frame$method)),
    height = 9,
    units = "in"
  )
}
