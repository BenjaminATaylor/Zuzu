#!/usr/bin/env Rscript

# Run one optional sample-size or sequencing-depth resample.
#
# Arguments use `--key value` pairs:
#   --mode sample|depth
#   --raw-counts PATH
#   --samplesheet PATH
#   --reflevel LABEL
#   --repeat INTEGER
#   --target INTEGER
#   --seed INTEGER
#   --enet true|false
#   --enet-alpha NUMBER
#   --enet-pairs INTEGER
#   --enet-pi NUMBER
#   --enet-folds INTEGER
#   --enet-nlambda INTEGER
#   --enet-lambda-min-ratio NUMBER
#   --enet-seed INTEGER
#   --output-prefix PATH
#
# Sample targets are samples per phenotype. Depth targets are counts per sample.
# Four files are written: summary, gene calls, sample allocations, and library
# metrics, each beginning with the requested output prefix.

usage = paste(
  "Usage: zuzu_resample.R",
  "--mode sample|depth",
  "--raw-counts counts.tsv",
  "--samplesheet samplesheet.csv",
  "--reflevel LABEL",
  "--repeat INTEGER",
  "--target INTEGER",
  "--seed INTEGER",
  "--enet true|false",
  "--enet-alpha NUMBER",
  "--enet-pairs INTEGER",
  "--enet-pi NUMBER",
  "--enet-folds INTEGER",
  "--enet-nlambda INTEGER",
  "--enet-lambda-min-ratio NUMBER",
  "--enet-seed INTEGER",
  "--output-prefix PATH"
)

# Parse the simple key-value command line
parse.arguments = function(arguments){
  if(length(arguments) == 1L && arguments[[1]] %in% c("--help", "-h")){
    cat(usage, "\n")
    quit(save = "no", status = 0L)
  }
  if(length(arguments) == 0L || length(arguments) %% 2L != 0L){
    stop(usage, call. = FALSE)
  }
  keys = arguments[seq.int(1L, length(arguments), by = 2L)]
  values = arguments[seq.int(2L, length(arguments), by = 2L)]
  if(any(!startsWith(keys, "--")) || anyDuplicated(keys)){
    stop(usage, call. = FALSE)
  }
  names(values) = substring(keys, 3L)
  as.list(values)
}

options = parse.arguments(commandArgs(trailingOnly = TRUE))
required.arguments = c(
  "mode", "raw-counts", "samplesheet", "reflevel",
  "repeat", "target", "seed", "enet", "enet-alpha", "enet-pairs",
  "enet-pi", "enet-folds", "enet-nlambda", "enet-lambda-min-ratio",
  "enet-seed", "output-prefix"
)
missing.arguments = setdiff(required.arguments, names(options))
unexpected.arguments = setdiff(names(options), required.arguments)
if(length(missing.arguments) > 0L || length(unexpected.arguments) > 0L){
  stop(sprintf(
    "Invalid arguments; missing: %s; unexpected: %s\n%s",
    if(length(missing.arguments) == 0L) "none" else
      paste(missing.arguments, collapse = ", "),
    if(length(unexpected.arguments) == 0L) "none" else
      paste(unexpected.arguments, collapse = ", "),
    usage
  ), call. = FALSE)
}

# Resolve and validate the requested settings
mode = options[["mode"]]
raw.counts.path = options[["raw-counts"]]
samplesheet.path = options[["samplesheet"]]
reflevel = options[["reflevel"]]
repeat.id = suppressWarnings(as.integer(options[["repeat"]]))
target = suppressWarnings(as.numeric(options[["target"]]))
data.seed = suppressWarnings(as.integer(options[["seed"]]))
enet.enabled = tolower(options[["enet"]])
enet.alpha = suppressWarnings(as.numeric(options[["enet-alpha"]]))
enet.pairs = suppressWarnings(as.integer(options[["enet-pairs"]]))
enet.pi = suppressWarnings(as.numeric(options[["enet-pi"]]))
enet.folds = suppressWarnings(as.integer(options[["enet-folds"]]))
enet.nlambda = suppressWarnings(as.integer(options[["enet-nlambda"]]))
enet.lambda.min.ratio = suppressWarnings(
  as.numeric(options[["enet-lambda-min-ratio"]])
)
enet.seed = suppressWarnings(as.integer(options[["enet-seed"]]))
output.prefix = options[["output-prefix"]]

if(!mode %in% c("sample", "depth")) stop("Mode must be sample or depth")
if(!enet.enabled %in% c("true", "false")) stop("--enet must be true or false")
enet.enabled = enet.enabled == "true"
if(any(!file.exists(c(raw.counts.path, samplesheet.path)))){
  stop("One or more input files do not exist")
}
if(is.na(repeat.id) || repeat.id < 1L) stop("Repeat must be a positive integer")
if(!is.finite(target) || target < 1 || target != floor(target)){
  stop("Target must be a positive integer")
}
if(is.na(data.seed) || data.seed < 0L) stop("Seed must be a non-negative integer")
if(is.na(enet.seed) || enet.seed < 0L){
  stop("Elastic-net seed must be a non-negative integer")
}
if(!is.finite(enet.alpha) || enet.alpha <= 0 || enet.alpha > 1){
  stop("Elastic-net alpha must be in (0, 1]")
}
if(is.na(enet.pairs) || enet.pairs < 1L){
  stop("Elastic-net pairs must be a positive integer")
}
if(!is.finite(enet.pi) || enet.pi <= 0.5 || enet.pi > 1){
  stop("Elastic-net pi must be in (0.5, 1]")
}
if(is.na(enet.folds) || enet.folds < 3L){
  stop("Elastic-net folds must be at least three")
}
if(is.na(enet.nlambda) || enet.nlambda < 2L){
  stop("Elastic-net nlambda must be at least two")
}
if(!is.finite(enet.lambda.min.ratio) ||
   enet.lambda.min.ratio <= 0 || enet.lambda.min.ratio >= 1){
  stop("Elastic-net lambda-min-ratio must be in (0, 1)")
}
target = as.numeric(target)
output.directory = dirname(output.prefix)
dir.create(output.directory, recursive = TRUE, showWarnings = FALSE)

# Read and align the sample information and raw count matrix
samplesheet = read.csv(
  samplesheet.path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if(!all(c("sample", "phenotype") %in% names(samplesheet))){
  stop("Sample sheet must contain sample and phenotype columns")
}
if(anyNA(samplesheet$sample) || anyDuplicated(samplesheet$sample)){
  stop("Sample sheet contains invalid sample identifiers")
}
samplesheet$phenotype = relevel(factor(samplesheet$phenotype), ref = reflevel)
if(nlevels(samplesheet$phenotype) != 2L){
  stop("Sample sheet phenotype must contain exactly two groups")
}
altlevel = setdiff(levels(samplesheet$phenotype), reflevel)
if(length(altlevel) != 1L) stop("Reference level is not valid for this dataset")

raw.frame = read.delim(
  raw.counts.path,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
missing.samples = setdiff(samplesheet$sample, names(raw.frame))
if(length(missing.samples) > 0L){
  stop(sprintf("Raw count matrix is missing %d samples", length(missing.samples)))
}
raw.counts = as.matrix(raw.frame[,samplesheet$sample,drop = FALSE])
storage.mode(raw.counts) = "double"
raw.counts = round(raw.counts)
if(any(!is.finite(raw.counts)) || any(raw.counts < 0)){
  stop("Raw count matrix contains missing, non-finite, or negative values")
}
if(anyDuplicated(row.names(raw.counts))){
  stop("Raw count matrix contains duplicate gene identifiers")
}
raw.genes = row.names(raw.counts)

# Draw an exact multivariate hypergeometric count vector
thin.count.vector = function(counts, requested.total){
  original.total = sum(counts)
  if(requested.total > original.total){
    stop(sprintf(
      "Requested depth %.0f exceeds a sample library of %.0f counts",
      requested.total, original.total
    ))
  }
  if(requested.total == original.total) return(counts)

  thinned = numeric(length(counts))
  remaining.draws = requested.total
  remaining.population = original.total
  if(length(counts) > 1L){
    for(index in seq_len(length(counts) - 1L)){
      this.count = counts[[index]]
      if(remaining.draws == 0 || this.count == 0){
        this.draw = 0
      } else if(remaining.draws == remaining.population){
        this.draw = this.count
      } else {
        this.draw = rhyper(
          1L,
          m = this.count,
          n = remaining.population - this.count,
          k = remaining.draws
        )
      }
      thinned[[index]] = this.draw
      remaining.draws = remaining.draws - this.draw
      remaining.population = remaining.population - this.count
    }
  }
  thinned[[length(counts)]] = remaining.draws
  if(sum(thinned) != requested.total || any(thinned < 0) || any(thinned > counts)){
    stop("Exact count thinning produced an invalid library")
  }
  thinned
}

# Generate the requested sample-size or depth-resampled dataset
set.seed(data.seed)
if(mode == "sample"){
  group.members = split(samplesheet$sample, samplesheet$phenotype)
  if(any(lengths(group.members) < target)){
    stop("Sample target exceeds the available size of at least one phenotype")
  }
  selected.ids = unlist(lapply(group.members, function(ids){
    sample(ids, size = target, replace = FALSE)
  }), use.names = FALSE)
  selected.ids = samplesheet$sample[samplesheet$sample %in% selected.ids]
  analysis.samples = samplesheet[
    match(selected.ids, samplesheet$sample),,
    drop = FALSE
  ]
  generated.counts = raw.counts[,selected.ids,drop = FALSE]
  original.library.sizes = colSums(generated.counts)
  prefilter.library.sizes = original.library.sizes
} else {
  analysis.samples = samplesheet
  generated.counts = apply(raw.counts, 2L, thin.count.vector, requested.total = target)
  if(is.null(dim(generated.counts))){
    generated.counts = matrix(
      generated.counts,
      ncol = 1L,
      dimnames = list(raw.genes, samplesheet$sample)
    )
  }
  row.names(generated.counts) = raw.genes
  colnames(generated.counts) = samplesheet$sample
  original.library.sizes = colSums(raw.counts)
  prefilter.library.sizes = colSums(generated.counts)
  if(any(prefilter.library.sizes != target)){
    stop("Depth thinning did not produce the exact requested sample totals")
  }
}

# Apply the common CPM and prevalence filter to the generated dataset
if(any(prefilter.library.sizes <= 0)){
  stop("Every generated sample must have a positive pre-filter library size")
}
minimum.samples = min(as.integer(table(analysis.samples$phenotype)))
at.least.one.cpm = sweep(
  generated.counts,
  MARGIN = 2L,
  STATS = prefilter.library.sizes / 1e6,
  FUN = ">="
)
keep.rows = rowSums(at.least.one.cpm) >= minimum.samples
if(!any(keep.rows)) stop("No genes pass the common CPM/prevalence filter")
analysis.counts = generated.counts[keep.rows,,drop = FALSE]
analysis.samples = analysis.samples[
  match(colnames(analysis.counts), analysis.samples$sample),,
  drop = FALSE
]
row.names(analysis.samples) = analysis.samples$sample
tested.genes = row.names(analysis.counts)

# Run DESeq2 using the basic pipeline model and contrast
run.deseq = function(counts, samples){
  if(!requireNamespace("DESeq2", quietly = TRUE)) stop("DESeq2 is not installed")
  dds = DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = samples,
    design = stats::as.formula("~phenotype")
  )
  dds = DESeq2::DESeq(dds, fitType = "parametric", betaPrior = FALSE, quiet = TRUE)
  result = DESeq2::results(
    dds,
    contrast = c("phenotype", altlevel, reflevel)
  )
  row.names(result)[!is.na(result$padj) & result$padj < 0.05]
}

# Run edgeR using the basic pipeline quasi-likelihood workflow
run.edger = function(counts, samples){
  if(!requireNamespace("edgeR", quietly = TRUE)) stop("edgeR is not installed")
  dge = edgeR::DGEList(counts = counts, group = samples$phenotype)
  dge = edgeR::normLibSizes(dge)
  design = stats::model.matrix(~phenotype, data = samples)
  phenotype.coef = which(attr(design, "assign") == 1L)
  contrast = rep(0, ncol(design))
  contrast[phenotype.coef] = 1
  dge = edgeR::estimateDisp(dge, design)
  fit = edgeR::glmQLFit(dge, design)
  test = edgeR::glmQLFTest(fit, contrast = contrast)
  adjusted = stats::p.adjust(test$table$PValue, method = "BH")
  row.names(test$table)[!is.na(adjusted) & adjusted < 0.05]
}

# Run per-gene Wilcoxon tests after TMM and CPM normalization
run.wilcoxon = function(counts, samples){
  if(!requireNamespace("edgeR", quietly = TRUE)) stop("edgeR is not installed")
  dge = edgeR::DGEList(counts = counts, group = samples$phenotype)
  dge = edgeR::calcNormFactors(dge, method = "TMM")
  normalized = edgeR::cpm(dge)
  reference.columns = samples$sample[samples$phenotype == reflevel]
  alternative.columns = samples$sample[samples$phenotype == altlevel]
  p.values = vapply(seq_len(nrow(normalized)), function(index){
    suppressWarnings(stats::wilcox.test(
      as.numeric(normalized[index,reference.columns]),
      as.numeric(normalized[index,alternative.columns])
    )$p.value)
  }, numeric(1))
  adjusted = stats::p.adjust(p.values, method = "BH")
  row.names(normalized)[!is.na(adjusted) & adjusted < 0.05]
}

deseq.calls = run.deseq(analysis.counts, analysis.samples)
edger.calls = run.edger(analysis.counts, analysis.samples)
wilcoxon.calls = run.wilcoxon(analysis.counts, analysis.samples)

# Locate and run the existing compact elastic-net selector
locate.enet.script = function(){
  full.arguments = commandArgs(trailingOnly = FALSE)
  file.argument = grep("^--file=", full.arguments, value = TRUE)
  if(length(file.argument) > 0L){
    this.script = sub("^--file=", "", file.argument[[1]])
    candidate = file.path(dirname(normalizePath(this.script)), "zuzu_enet.R")
    if(file.exists(candidate)) return(candidate)
  }
  candidate = Sys.which("zuzu_enet.R")
  if(nzchar(candidate)) return(unname(candidate))
  stop("Could not locate bin/zuzu_enet.R")
}

run.enet = function(counts, samples, library.sizes){
  enet.script = locate.enet.script()
  work.directory = tempfile(pattern = "zuzu_enet_")
  dir.create(work.directory)
  on.exit(unlink(work.directory, recursive = TRUE), add = TRUE)

  write.csv(samples, file.path(work.directory, "samplesheet.csv"), row.names = FALSE)
  write.csv(counts, file.path(work.directory, "counts.csv"))
  write.csv(
    data.frame(
      sample = colnames(counts),
      library_size = as.numeric(library.sizes[colnames(counts)])
    ),
    file.path(work.directory, "library_sizes.csv"),
    row.names = FALSE
  )

  command.arguments = c(
    shQuote(enet.script),
    "--mode", "compact",
    "--samplesheet", "samplesheet.csv",
    "--counts", "counts.csv",
    "--reflevel", shQuote(reflevel),
    "--library-sizes", "library_sizes.csv",
    "--alpha", as.character(enet.alpha),
    "--pairs", as.character(enet.pairs),
    "--pi", as.character(enet.pi),
    "--folds", as.character(enet.folds),
    "--nlambda", as.character(enet.nlambda),
    "--lambda-min-ratio", as.character(enet.lambda.min.ratio),
    "--seed", as.character(enet.seed)
  )

  previous.directory = getwd()
  status = tryCatch({
    setwd(work.directory)
    system2("Rscript", command.arguments)
  }, finally = setwd(previous.directory))
  if(!identical(status, 0L)) stop("Elastic-net selector failed")

  stable.path = file.path(work.directory, "enet_stable_features.csv")
  summary.path = file.path(work.directory, "enet_summary.csv")
  if(!file.exists(stable.path) || !file.exists(summary.path)){
    stop("Elastic-net selector did not produce its compact outputs")
  }
  stable = read.csv(stable.path, check.names = FALSE, stringsAsFactors = FALSE)
  summary = read.csv(summary.path, check.names = FALSE, stringsAsFactors = FALSE)
  if(!"gene" %in% names(stable)){
    stop("Elastic-net stable-feature output does not contain a gene column")
  }
  list(genes = as.character(stable$gene), summary = summary)
}

if(enet.enabled){
  enet.result = run.enet(
    analysis.counts,
    analysis.samples,
    prefilter.library.sizes
  )
  enet.calls = enet.result$genes
  enet.status = "ok"
  enet.maximum.probability = if(
    "maximum_selection_probability" %in% names(enet.result$summary)
  ) enet.result$summary$maximum_selection_probability[[1]] else NA_real_
} else {
  enet.calls = character()
  enet.status = "disabled"
  enet.maximum.probability = NA_real_
}

# Assemble one complete gene-selection table over the raw input universe
elastic.calls = if(enet.enabled){
  raw.genes %in% enet.calls
} else {
  rep(NA, length(raw.genes))
}
gene.results = data.frame(
  gene = raw.genes,
  DESeq2 = raw.genes %in% deseq.calls,
  edgeR = raw.genes %in% edger.calls,
  Wilcoxon = raw.genes %in% wilcoxon.calls,
  `Elastic net` = elastic.calls,
  tested = raw.genes %in% tested.genes,
  repeat_id = repeat.id,
  target = target,
  resampling_type = mode,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
names(gene.results)[names(gene.results) == "repeat_id"] = "repeat"

# Record the selected samples used by sample-size resampling
if(mode == "sample"){
  allocations = data.frame(
    resampling_type = mode,
    repeat_id = repeat.id,
    target = target,
    sample = analysis.samples$sample,
    phenotype = as.character(analysis.samples$phenotype),
    stringsAsFactors = FALSE
  )
} else {
  allocations = data.frame(
    resampling_type = character(),
    repeat_id = integer(),
    target = numeric(),
    sample = character(),
    phenotype = character(),
    stringsAsFactors = FALSE
  )
}
names(allocations)[names(allocations) == "repeat_id"] = "repeat"

# Record library depths before thinning, before filtering, and after filtering
library.metrics = data.frame(
  resampling_type = mode,
  repeat_id = repeat.id,
  target = target,
  sample = analysis.samples$sample,
  phenotype = as.character(analysis.samples$phenotype),
  original_depth = as.numeric(original.library.sizes[analysis.samples$sample]),
  target_depth = if(mode == "depth") rep(target, nrow(analysis.samples)) else
    rep(NA_real_, nrow(analysis.samples)),
  thinned_prefilter_depth = as.numeric(
    prefilter.library.sizes[analysis.samples$sample]
  ),
  postfilter_depth = as.numeric(colSums(analysis.counts)),
  stringsAsFactors = FALSE
)
names(library.metrics)[names(library.metrics) == "repeat_id"] = "repeat"

# Record one compact dataset-level summary row
group.sizes = table(analysis.samples$phenotype)
summary.result = data.frame(
  resampling_type = mode,
  repeat_id = repeat.id,
  target = target,
  target_definition = if(mode == "sample") "samples_per_phenotype" else
    "counts_per_sample",
  generation_seed = data.seed,
  samples_per_group = if(mode == "sample") target else NA_real_,
  target_depth = if(mode == "depth") target else NA_real_,
  reference = reflevel,
  alternative = altlevel,
  reference_n = as.integer(group.sizes[[reflevel]]),
  alternative_n = as.integer(group.sizes[[altlevel]]),
  total_n = nrow(analysis.samples),
  n_input_genes = length(raw.genes),
  n_tested = length(tested.genes),
  filtered_genes = length(raw.genes) - length(tested.genes),
  minimum_samples_for_one_cpm = minimum.samples,
  DESeq2 = length(deseq.calls),
  edgeR = length(edger.calls),
  Wilcoxon = length(wilcoxon.calls),
  `Elastic net` = if(enet.enabled) length(enet.calls) else NA_integer_,
  enet_enabled = enet.enabled,
  enet_status = enet.status,
  enet_alpha = enet.alpha,
  enet_pairs = enet.pairs,
  enet_pi = enet.pi,
  enet_folds = enet.folds,
  enet_nlambda = enet.nlambda,
  enet_lambda_min_ratio = enet.lambda.min.ratio,
  enet_seed = enet.seed,
  enet_maximum_selection_probability = enet.maximum.probability,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
names(summary.result)[names(summary.result) == "repeat_id"] = "repeat"

# Write the four requested resampling artifacts
write.csv(
  summary.result,
  paste0(output.prefix, "_summary.csv"),
  row.names = FALSE,
  na = "NA"
)
write.csv(
  gene.results,
  paste0(output.prefix, "_genes.csv"),
  row.names = FALSE,
  na = "NA"
)
write.csv(
  allocations,
  paste0(output.prefix, "_allocations.csv"),
  row.names = FALSE,
  na = "NA"
)
write.csv(
  library.metrics,
  paste0(output.prefix, "_library_metrics.csv"),
  row.names = FALSE,
  na = "NA"
)
