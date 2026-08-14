#!/usr/bin/env Rscript

# Keep numerical libraries single-threaded so that Nextflow can parallelise datasets
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)

# Check that the small set of runtime dependencies is available
required.packages = c("glmnet", "ggplot2", "jsonlite")
missing.packages = required.packages[!vapply(
  required.packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)]
if(length(missing.packages) > 0L){
  stop(
    "Missing required R package(s): ",
    paste(missing.packages, collapse = ", "),
    call. = FALSE
  )
}

# Parse the deliberately simple named-pair command line interface
parse.arguments = function(arguments){
  expected = c(
    "mode", "samplesheet", "counts", "reflevel", "library-sizes",
    "alpha", "pairs", "pi", "folds", "nlambda", "lambda-min-ratio",
    "seed"
  )
  if(length(arguments) == 0L || length(arguments) %% 2L != 0L){
    stop(
      "Arguments must be supplied as named pairs, for example --mode basic",
      call. = FALSE
    )
  }
  keys = arguments[seq.int(1L, length(arguments), by = 2L)]
  values = arguments[seq.int(2L, length(arguments), by = 2L)]
  if(any(!startsWith(keys, "--"))){
    stop("Every argument name must begin with --", call. = FALSE)
  }
  keys = substring(keys, 3L)
  if(any(keys == "")) stop("An empty argument name was supplied", call. = FALSE)
  if(anyDuplicated(keys)){
    stop(
      "Duplicate argument(s): ",
      paste(unique(keys[duplicated(keys)]), collapse = ", "),
      call. = FALSE
    )
  }
  unknown = setdiff(keys, expected)
  missing = setdiff(expected, keys)
  if(length(unknown) > 0L){
    stop("Unknown argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  if(length(missing) > 0L){
    stop("Missing argument(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  stats::setNames(as.list(values), keys)
}

# Convert and validate scalar numerical settings without silent truncation
parse.integer = function(value, label, minimum = 1L){
  parsed = suppressWarnings(as.numeric(value))
  if(length(parsed) != 1L || !is.finite(parsed) || parsed != floor(parsed) ||
     parsed < minimum || parsed > .Machine$integer.max){
    stop(label, " must be an integer >= ", minimum, call. = FALSE)
  }
  as.integer(parsed)
}

parse.number = function(value, label){
  parsed = suppressWarnings(as.numeric(value))
  if(length(parsed) != 1L || !is.finite(parsed)){
    stop(label, " must be a finite number", call. = FALSE)
  }
  parsed
}

arguments = parse.arguments(commandArgs(trailingOnly = TRUE))
mode = arguments[["mode"]]
if(!mode %in% c("basic", "compact")){
  stop("--mode must be basic or compact", call. = FALSE)
}

samplesheet.path = arguments[["samplesheet"]]
counts.path = arguments[["counts"]]
reflevel = arguments[["reflevel"]]
library.sizes.argument = arguments[["library-sizes"]]
alpha = parse.number(arguments[["alpha"]], "--alpha")
n.pairs = parse.integer(arguments[["pairs"]], "--pairs")
pi.threshold = parse.number(arguments[["pi"]], "--pi")
requested.folds = parse.integer(arguments[["folds"]], "--folds", minimum = 3L)
nlambda = parse.integer(arguments[["nlambda"]], "--nlambda", minimum = 3L)
lambda.min.ratio = parse.number(
  arguments[["lambda-min-ratio"]],
  "--lambda-min-ratio"
)
master.seed = parse.integer(arguments[["seed"]], "--seed", minimum = 0L)

if(n.pairs > .Machine$integer.max %/% 2L){
  stop("--pairs is too large", call. = FALSE)
}
if(alpha <= 0 || alpha > 1){
  stop("--alpha must be greater than 0 and no greater than 1", call. = FALSE)
}
if(pi.threshold <= 0.5 || pi.threshold > 1){
  stop("--pi must be greater than 0.5 and no greater than 1", call. = FALSE)
}
if(lambda.min.ratio <= 0 || lambda.min.ratio >= 1){
  stop("--lambda-min-ratio must be greater than 0 and less than 1", call. = FALSE)
}
if(!nzchar(reflevel)) stop("--reflevel must not be empty", call. = FALSE)
if(!file.exists(samplesheet.path)) stop("Sample sheet does not exist", call. = FALSE)
if(!file.exists(counts.path)) stop("Count matrix does not exist", call. = FALSE)
if(library.sizes.argument != "NONE" && !file.exists(library.sizes.argument)){
  stop("Library-size sidecar does not exist", call. = FALSE)
}

# Read comma- or tab-delimited inputs without changing sample identifiers
read.input.table = function(path){
  if(grepl("[.](tsv|txt)([.]gz)?$", path, ignore.case = TRUE)){
    utils::read.delim(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      comment.char = "",
      quote = "\""
    )
  } else {
    utils::read.csv(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      comment.char = ""
    )
  }
}

# Load and align the phenotype and count matrices
samples = read.input.table(samplesheet.path)
if(!all(c("sample", "phenotype") %in% names(samples))){
  stop("Sample sheet must contain sample and phenotype columns", call. = FALSE)
}
if(nrow(samples) == 0L) stop("Sample sheet contains no samples", call. = FALSE)
if(anyNA(samples$sample) || any(samples$sample == "") || anyDuplicated(samples$sample)){
  stop("Sample identifiers must be present and unique", call. = FALSE)
}
if(anyNA(samples$phenotype) || any(samples$phenotype == "")){
  stop("Every sample must have a phenotype", call. = FALSE)
}
samples$sample = as.character(samples$sample)
samples$phenotype = as.character(samples$phenotype)
phenotype.levels = unique(samples$phenotype)
if(length(phenotype.levels) != 2L || !reflevel %in% phenotype.levels){
  stop(
    "The sample sheet must contain exactly two phenotypes including --reflevel",
    call. = FALSE
  )
}
altlevel = setdiff(phenotype.levels, reflevel)
if(length(altlevel) != 1L) stop("Could not resolve the alternative phenotype", call. = FALSE)
sampling.levels = sort(phenotype.levels, method = "radix")

counts.frame = read.input.table(counts.path)
if(nrow(counts.frame) == 0L || ncol(counts.frame) < 2L){
  stop("Count matrix must contain genes and sample columns", call. = FALSE)
}
if(anyDuplicated(names(counts.frame))){
  stop("Count matrix contains duplicate column names", call. = FALSE)
}
gene.column = names(counts.frame)[[1L]]
gene.ids = as.character(counts.frame[[1L]])
if(anyNA(gene.ids) || any(gene.ids == "") || anyDuplicated(gene.ids)){
  stop("Gene identifiers must be present and unique", call. = FALSE)
}
missing.samples = setdiff(samples$sample, names(counts.frame))
allowed.annotations = c(gene.column, "gene_name")
unexpected.columns = setdiff(names(counts.frame), c(allowed.annotations, samples$sample))
if(length(missing.samples) > 0L || length(unexpected.columns) > 0L){
  stop(
    paste0(
      "Count/sample columns do not align; missing samples: ",
      if(length(missing.samples) == 0L) "none" else
        paste(missing.samples, collapse = ", "),
      "; unexpected columns: ",
      if(length(unexpected.columns) == 0L) "none" else
        paste(unexpected.columns, collapse = ", ")
    ),
    call. = FALSE
  )
}
if(any(!vapply(counts.frame[,samples$sample,drop = FALSE], is.numeric, logical(1)))){
  stop("Every count sample column must be numeric", call. = FALSE)
}

counts = as.matrix(counts.frame[,samples$sample,drop = FALSE])
storage.mode(counts) = "double"
row.names(counts) = gene.ids
if(any(!is.finite(counts)) || any(counts < 0)){
  stop("Count matrix contains missing, non-finite, or negative values", call. = FALSE)
}
counts = round(counts)
gene.names = if("gene_name" %in% names(counts.frame)){
  as.character(counts.frame$gene_name)
} else {
  gene.ids
}
gene.names[is.na(gene.names) | gene.names == ""] = gene.ids[
  is.na(gene.names) | gene.names == ""
]

# Resolve library sizes from the pre-filter sidecar or the supplied matrix
if(library.sizes.argument == "NONE"){
  library.sizes = colSums(counts)
  library.size.provenance = "column sums over every gene in the supplied matrix"
  library.sizes.path = NULL
} else {
  library.frame = read.input.table(library.sizes.argument)
  if(!all(c("sample", "library_size") %in% names(library.frame))){
    stop(
      "Library-size sidecar must contain sample and library_size columns",
      call. = FALSE
    )
  }
  if(anyNA(library.frame$sample) || anyDuplicated(library.frame$sample) ||
     !setequal(as.character(library.frame$sample), samples$sample)){
    stop("Library-size samples do not align with the sample sheet", call. = FALSE)
  }
  if(!is.numeric(library.frame$library_size)){
    stop("Library sizes must be numeric", call. = FALSE)
  }
  library.frame = library.frame[
    match(samples$sample, as.character(library.frame$sample)),,
    drop = FALSE
  ]
  library.sizes = as.numeric(library.frame$library_size)
  names(library.sizes) = samples$sample
  library.size.provenance = "supplied pre-filter library-size sidecar"
  library.sizes.path = normalizePath(library.sizes.argument, mustWork = TRUE)
}
if(any(!is.finite(library.sizes)) || any(library.sizes <= 0)){
  stop("Every library size must be finite and positive", call. = FALSE)
}
if(any(library.sizes + sqrt(.Machine$double.eps) < colSums(counts))){
  stop(
    "A supplied library size is smaller than its retained count-column sum",
    call. = FALSE
  )
}

# Transform counts to log2(CPM + 1) and encode the requested contrast
cpm = sweep(counts, 2L, library.sizes / 1e6, FUN = "/")
x = t(log2(cpm + 1))
colnames(x) = gene.ids
rownames(x) = samples$sample
y = as.integer(samples$phenotype == altlevel)
if(!identical(rownames(x), samples$sample)){
  stop("Internal sample alignment failed", call. = FALSE)
}
n.genes = ncol(x)
if(n.genes == 0L) stop("The analysis universe contains no genes", call. = FALSE)
class.counts.full = table(factor(y, levels = c(0L, 1L)))
if(any(class.counts.full == 0L)) stop("Both phenotypes require samples", call. = FALSE)

# Generate complementary phenotype-stratified half-samples before fitting
set.seed(master.seed)
class.indices = split(
  seq_along(y),
  factor(samples$phenotype, levels = sampling.levels)
)
total.fits = 2L * n.pairs
design = vector("list", total.fits)
fit.index = 0L
for(pair.id in seq_len(n.pairs)){
  half.a = integer()
  for(class.position in seq_along(class.indices)){
    indices = class.indices[[class.position]]
    smaller.size = floor(length(indices) / 2)
    larger.size = ceiling(length(indices) / 2)
    use.larger = length(indices) %% 2L == 1L &&
      (pair.id + class.position) %% 2L == 1L
    draw.size = if(use.larger) larger.size else smaller.size
    if(draw.size == 0L){
      stop("A phenotype is too small to form complementary halves", call. = FALSE)
    }
    half.a = c(
      half.a,
      sample(indices, size = draw.size, replace = FALSE)
    )
  }
  halves = list(A = sort(half.a), B = setdiff(seq_along(y), half.a))
  if(length(intersect(halves$A, halves$B)) != 0L ||
     !identical(sort(c(halves$A, halves$B)), seq_along(y))){
    stop("Internal complementary-half construction failed", call. = FALSE)
  }
  for(half.id in names(halves)){
    fit.index = fit.index + 1L
    design[[fit.index]] = list(
      fit = fit.index,
      pair = pair.id,
      half = half.id,
      indices = halves[[half.id]]
    )
  }
}

# Resolve one safe fold count shared by all half-sample fits
smallest.half.class = min(vapply(
  design,
  function(item){
    min(table(factor(y[item$indices], levels = c(0L, 1L))))
  },
  integer(1)
))
smallest.grouped.cv.limit = min(vapply(
  design,
  function(item) as.integer(floor(length(item$indices) / 3L)),
  integer(1)
))
resolved.folds = min(
  requested.folds,
  smallest.half.class,
  smallest.grouped.cv.limit
)
if(resolved.folds < 3L){
  stop(
    paste0(
      "At least three safe inner folds are required after half-sampling; ",
      "increase the number of samples in the smaller phenotype"
    ),
    call. = FALSE
  )
}

# Assign deterministic phenotype-stratified folds within every half-sample
for(index in seq_along(design)){
  indices = design[[index]]$indices
  fold.id = integer(length(indices))
  for(phenotype.level in sampling.levels){
    local = which(samples$phenotype[indices] == phenotype.level)
    fold.id[local] = sample(
      rep(seq_len(resolved.folds), length.out = length(local))
    )
  }
  if(!all(vapply(
    seq_len(resolved.folds),
    function(fold){
      length(unique(y[indices][fold.id == fold])) == 2L
    },
    logical(1)
  ))){
    stop("Internal stratified-fold construction failed", call. = FALSE)
  }
  design[[index]]$fold.id = fold.id
}

# Fit one uncapped cross-validated elastic-net selector to a half-sample
fit.base.selector = function(x.half, y.half, fold.id){
  class.counts = table(factor(y.half, levels = c(0L, 1L)))
  if(any(class.counts == 0L)) stop("A half-sample contains only one phenotype")
  variances = apply(x.half, 2L, stats::var)
  variable = is.finite(variances) & variances > 0
  if(!any(variable)) stop("A half-sample contains no variable genes")
  fit.x = x.half[,variable,drop = FALSE]
  weights = length(y.half) /
    (2 * as.numeric(class.counts[as.character(y.half)]))

  cv.fit = glmnet::cv.glmnet(
    x = fit.x,
    y = y.half,
    family = "binomial",
    alpha = alpha,
    weights = weights,
    foldid = fold.id,
    type.measure = "deviance",
    alignment = "lambda",
    grouped = TRUE,
    standardize = TRUE,
    intercept = TRUE,
    nlambda = nlambda,
    lambda.min.ratio = lambda.min.ratio,
    parallel = FALSE,
    keep = FALSE
  )
  coefficients = as.matrix(stats::coef(cv.fit, s = "lambda.1se"))[-1L,1L]
  selected = rep(FALSE, ncol(x.half))
  selected[which(variable)[which(coefficients != 0)]] = TRUE
  path.minimum = min(cv.fit$lambda)

  list(
    selected = selected,
    selected.count = sum(selected),
    lambda.1se = unname(cv.fit$lambda.1se),
    lambda.min = unname(cv.fit$lambda.min),
    lambda.max = unname(max(cv.fit$lambda)),
    path.lambda.min = unname(path.minimum),
    lambda.1se.at.path.min =
      abs(cv.fit$lambda.1se - path.minimum) <=
        sqrt(.Machine$double.eps) * max(1, abs(path.minimum)),
    path.length = length(cv.fit$lambda),
    zero.variance.count = sum(!variable),
    class.0 = unname(class.counts[[1L]]),
    class.1 = unname(class.counts[[2L]])
  )
}

# Run all complementary-half fits and retain complete diagnostics
started = proc.time()[["elapsed"]]
selection.matrix = matrix(
  FALSE,
  nrow = total.fits,
  ncol = n.genes,
  dimnames = list(as.character(seq_len(total.fits)), gene.ids)
)
diagnostic.rows = vector("list", total.fits)
allocation.rows = if(mode == "basic") vector("list", total.fits) else NULL

for(index in seq_along(design)){
  fit.design = design[[index]]
  indices = fit.design$indices
  result = tryCatch(
    fit.base.selector(
      x[indices,,drop = FALSE],
      y[indices],
      fit.design$fold.id
    ),
    error = function(error){
      stop(
        sprintf(
          "Elastic-net fit failed at pair %d half %s: %s",
          fit.design$pair,
          fit.design$half,
          conditionMessage(error)
        ),
        call. = FALSE
      )
    }
  )
  selection.matrix[index,] = result$selected
  diagnostic.rows[[index]] = data.frame(
    fit = fit.design$fit,
    pair = fit.design$pair,
    half = fit.design$half,
    n_samples = length(indices),
    n_reference = result$class.0,
    n_alternative = result$class.1,
    requested_folds = requested.folds,
    resolved_folds = resolved.folds,
    selected_count = result$selected.count,
    lambda_1se = result$lambda.1se,
    lambda_min = result$lambda.min,
    lambda_max = result$lambda.max,
    path_lambda_min = result$path.lambda.min,
    lambda_1se_ratio = result$lambda.1se / result$lambda.max,
    lambda_1se_at_path_min = result$lambda.1se.at.path.min,
    path_length = result$path.length,
    zero_variance_count = result$zero.variance.count,
    status = "ok",
    stringsAsFactors = FALSE
  )
  if(mode == "basic"){
    allocation.rows[[index]] = data.frame(
      fit = fit.design$fit,
      pair = fit.design$pair,
      half = fit.design$half,
      sample = samples$sample[indices],
      phenotype = samples$phenotype[indices],
      class_value = y[indices],
      inner_fold = fit.design$fold.id,
      stringsAsFactors = FALSE
    )
  }
}
selector.elapsed = proc.time()[["elapsed"]] - started
fit.diagnostics = do.call(rbind, diagnostic.rows)

# Convert selection counts into one deterministic complete gene ranking
selection.counts = as.integer(colSums(selection.matrix))
selection.probabilities = selection.counts / total.fits
threshold.count = as.integer(ceiling(pi.threshold * total.fits - 1e-12))
stable = selection.counts >= threshold.count
ordering = order(-selection.probabilities, gene.ids, method = "radix")
gene.table = data.frame(
  gene = gene.ids[ordering],
  gene_name = gene.names[ordering],
  selection_count = selection.counts[ordering],
  total_fits = total.fits,
  selection_probability = selection.probabilities[ordering],
  stable_feature = stable[ordering],
  rank = seq_along(ordering),
  stringsAsFactors = FALSE
)
stable.table = gene.table[gene.table$stable_feature,,drop = FALSE]
status = if(nrow(stable.table) == 0L) "ok_empty" else "ok"

# Write the compact outputs required for every selector run
summary.table = data.frame(
  method = "Stability-selected elastic net",
  status = status,
  n_tested = n.genes,
  n_called = nrow(stable.table),
  stable_feature_count = nrow(stable.table),
  maximum_selection_probability = max(selection.probabilities),
  mean_base_selected = mean(fit.diagnostics$selected_count),
  median_base_selected = stats::median(fit.diagnostics$selected_count),
  minimum_base_selected = min(fit.diagnostics$selected_count),
  maximum_base_selected = max(fit.diagnostics$selected_count),
  empty_base_fits = sum(fit.diagnostics$selected_count == 0L),
  alpha = alpha,
  complementary_pairs = n.pairs,
  total_half_sample_fits = total.fits,
  pi_threshold = pi.threshold,
  threshold_count = threshold.count,
  requested_folds = requested.folds,
  resolved_folds = resolved.folds,
  nlambda = nlambda,
  lambda_min_ratio = lambda.min.ratio,
  seed = master.seed,
  q_status = "unset_and_uncapped",
  pfer_claim = "none",
  elapsed_seconds = unname(selector.elapsed),
  stringsAsFactors = FALSE
)

utils::write.csv(
  stable.table,
  "enet_stable_features.csv",
  row.names = FALSE
)
utils::write.csv(
  summary.table,
  "enet_summary.csv",
  row.names = FALSE
)
utils::write.csv(
  fit.diagnostics,
  "enet_fit_diagnostics.csv",
  row.names = FALSE
)

# Retain complete provenance without making a PFER or gene-wise error-control claim
metadata = list(
  status = status,
  mode = mode,
  interpretation = paste0(
    "Stable features are genes repeatedly selected by an uncapped sparse additive ",
    "logistic model under complementary half-sample perturbations. Selection ",
    "probabilities are resampling frequencies, not p-values or posterior probabilities."
  ),
  inputs = list(
    samplesheet = normalizePath(samplesheet.path, mustWork = TRUE),
    counts = normalizePath(counts.path, mustWork = TRUE),
    library_sizes = library.sizes.path,
    samplesheet_md5 = unname(tools::md5sum(samplesheet.path)),
    counts_md5 = unname(tools::md5sum(counts.path)),
    library_sizes_md5 = if(is.null(library.sizes.path)) NULL else
      unname(tools::md5sum(library.sizes.argument)),
    samples = nrow(samples),
    genes = n.genes,
    reference = reflevel,
    alternative = altlevel,
    reference_samples = unname(class.counts.full[[1L]]),
    alternative_samples = unname(class.counts.full[[2L]])
  ),
  preprocessing = list(
    counts_rounded = TRUE,
    additional_gene_filter = NULL,
    library_size_provenance = library.size.provenance,
    transformation = "log2(CPM + 1)",
    predictor_standardization = "within each half-sample fit"
  ),
  selector = list(
    family = "binomial",
    alpha = alpha,
    complementary_pairs = n.pairs,
    total_half_sample_fits = total.fits,
    sampling_phenotype_order = sampling.levels,
    contrast_invariant_randomization = TRUE,
    q = NULL,
    q_status = "unset_and_uncapped",
    base_model_rule = "within-half stratified CV at lambda.1se",
    requested_cv_folds = requested.folds,
    resolved_cv_folds = resolved.folds,
    fold_resolution_rule = paste0(
      "min(requested folds, smallest half-class count, ",
      "smallest floor(half-sample size / 3)); require at least 3"
    ),
    cv_measure = "binomial deviance",
    alignment = "lambda",
    grouped = TRUE,
    standardize = TRUE,
    intercept = TRUE,
    class_weights = "n_half / (2 * n_class_in_half)",
    nlambda_requested = nlambda,
    lambda_min_ratio = lambda.min.ratio,
    lambda_rule = "lambda.1se",
    pi_threshold = pi.threshold,
    threshold_count = threshold.count,
    threshold_inclusive = TRUE,
    pfer_bound = NULL,
    pfer_claim = "none",
    error_control_caution = paste0(
      "No q cap or nominal PFER bound is used. The stability threshold is not ",
      "gene-wise FDR control."
    ),
    seed = master.seed
  ),
  results = list(
    stable_feature_count = nrow(stable.table),
    maximum_selection_probability = max(selection.probabilities),
    mean_base_selected = mean(fit.diagnostics$selected_count),
    minimum_base_selected = min(fit.diagnostics$selected_count),
    maximum_base_selected = max(fit.diagnostics$selected_count),
    empty_base_fits = sum(fit.diagnostics$selected_count == 0L)
  ),
  execution = list(
    elapsed_seconds = unname(selector.elapsed),
    glmnet_parallel = FALSE,
    numerical_library_threads = 1L
  ),
  software = list(
    R = R.version.string,
    glmnet = as.character(utils::packageVersion("glmnet")),
    ggplot2 = as.character(utils::packageVersion("ggplot2")),
    jsonlite = as.character(utils::packageVersion("jsonlite"))
  )
)
jsonlite::write_json(
  metadata,
  "enet_metadata.json",
  pretty = TRUE,
  auto_unbox = TRUE,
  null = "null"
)

# Retain the full result and inspectable diagnostics for a baseline run
if(mode == "basic"){
  utils::write.csv(
    gene.table,
    "enet_gene_stability.csv",
    row.names = FALSE
  )
  allocations = do.call(rbind, allocation.rows)
  utils::write.csv(
    allocations,
    "enet_allocations_and_inner_folds.csv",
    row.names = FALSE
  )

  top.count = min(40L, nrow(gene.table))
  top.table = gene.table[seq_len(top.count),,drop = FALSE]
  top.table$gene_label = factor(
    top.table$gene,
    levels = rev(top.table$gene)
  )
  top.plot = ggplot2::ggplot(
    top.table,
    ggplot2::aes(
      x = gene_label,
      y = selection_probability,
      fill = stable_feature
    )
  ) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_hline(
      yintercept = pi.threshold,
      linetype = "dashed",
      colour = "#B2182B"
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = c(`FALSE` = "grey70", `TRUE` = "#0072B2"),
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.1),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      title = "Most frequently selected elastic-net features",
      subtitle = sprintf(
        "%d stable features at pi >= %.2f",
        nrow(stable.table),
        pi.threshold
      ),
      x = "Gene",
      y = "Selection probability",
      fill = "Stable feature"
    ) +
    ggplot2::theme_bw()
  ggplot2::ggsave(
    "enet_top_features.png",
    top.plot,
    width = 10,
    height = max(6, 0.18 * top.count + 2),
    units = "in",
    dpi = 160,
    bg = "white"
  )

  base.size.plot = ggplot2::ggplot(
    fit.diagnostics,
    ggplot2::aes(x = selected_count)
  ) +
    ggplot2::geom_histogram(
      binwidth = 1,
      boundary = -0.5,
      colour = "grey20",
      fill = "#56B4E9"
    ) +
    ggplot2::geom_vline(
      xintercept = mean(fit.diagnostics$selected_count),
      linetype = "dashed",
      colour = "#B2182B"
    ) +
    ggplot2::labs(
      title = "Elastic-net base-model sizes",
      subtitle = sprintf(
        "%d half-sample fits; mean %.1f selected genes",
        total.fits,
        mean(fit.diagnostics$selected_count)
      ),
      x = "Genes with non-zero coefficients",
      y = "Half-sample fits"
    ) +
    ggplot2::theme_bw()
  ggplot2::ggsave(
    "enet_base_model_sizes.png",
    base.size.plot,
    width = 9,
    height = 6,
    units = "in",
    dpi = 160,
    bg = "white"
  )

  top.genes = top.table$gene
  matrix.long = as.data.frame(
    as.table(selection.matrix[,top.genes,drop = FALSE]),
    stringsAsFactors = FALSE
  )
  names(matrix.long) = c("fit", "gene", "selected")
  matrix.long$fit = as.integer(as.character(matrix.long$fit))
  matrix.long$selected = as.logical(matrix.long$selected)
  matrix.long$gene = factor(matrix.long$gene, levels = rev(top.genes))
  selection.matrix.plot = ggplot2::ggplot(
    matrix.long,
    ggplot2::aes(x = fit, y = gene, fill = selected)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(
      values = c(`FALSE` = "white", `TRUE` = "#0072B2"),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "Top-feature selections across complementary halves",
      x = "Half-sample fit",
      y = "Gene",
      fill = "Selected"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  ggplot2::ggsave(
    "enet_top_feature_selection_matrix.png",
    selection.matrix.plot,
    width = 12,
    height = max(6, 0.18 * top.count + 2),
    units = "in",
    dpi = 160,
    bg = "white"
  )
}

message(sprintf(
  paste0(
    "Elastic net complete: %d/%d stable features at pi %.3f; ",
    "%d half-sample fits with %d inner folds"
  ),
  nrow(stable.table),
  n.genes,
  pi.threshold,
  total.fits,
  resolved.folds
))
