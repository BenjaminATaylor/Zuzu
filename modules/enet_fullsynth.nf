process ENET_FULLSYNTH {

  label 'mid_job'
  tag "replicate=${x},ref=${refnum},alt=${altnum},depth=${depth}"

  publishDir "$params.outdir/enet", pattern: 'enet_fullsynth_*'

  input:
  tuple val(depth),
        path(truedegs),
        path(samplesheet),
        path(countsframe),
        val(refnum),
        val(altnum),
        val(x)

  output:
  path "outframe.csv", emit: result
  tuple val(depth),
        val(refnum),
        val(altnum),
        val(x),
        path("enet_fullsynth_*_stable_features.csv"),
        path("enet_fullsynth_*_summary.csv"),
        path("enet_fullsynth_*_fit_diagnostics.csv"),
        path("enet_fullsynth_*_metadata.json"),
        emit: diagnostics

  script:
  depth_token = depth.toLong()
  task_seed = params.enet_seed.toInteger() + x.toInteger()
  output_prefix = "enet_fullsynth_depth_${depth_token}_ref_${refnum}_alt_${altnum}_replicate_${x}"
  """
  # Run the complete selector without exposing the synthetic truth set
  zuzu_enet.R \
    --mode compact \
    --samplesheet "$samplesheet" \
    --counts "$countsframe" \
    --reflevel '${params.reflevel}' \
    --library-sizes NONE \
    --alpha '${params.enet_alpha}' \
    --pairs '${params.enet_pairs}' \
    --pi '${params.enet_pi}' \
    --folds '${params.enet_folds}' \
    --nlambda '${params.enet_nlambda}' \
    --lambda-min-ratio '${params.enet_lambda_min_ratio}' \
    --seed '${task_seed}'

  # Give every retained synthetic artifact a unique name
  mv enet_stable_features.csv \
    "${output_prefix}_stable_features.csv"
  mv enet_summary.csv \
    "${output_prefix}_summary.csv"
  mv enet_fit_diagnostics.csv \
    "${output_prefix}_fit_diagnostics.csv"
  mv enet_metadata.json \
    "${output_prefix}_metadata.json"

  # Compare stable features with the simulator truth over the complete universe
  Rscript - <<'RSCRIPT'
  countsframe = read.csv(
    "$countsframe",
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  stable.frame = read.csv(
    "${output_prefix}_stable_features.csv",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  truedegs = get(load("$truedegs"))

  # Validate the fixed synthetic evaluation universe
  allgenes = row.names(countsframe)
  if(anyDuplicated(allgenes)){
    stop("Synthetic count matrix contains duplicate gene identifiers")
  }
  if(anyDuplicated(truedegs)){
    stop("Synthetic truth set contains duplicate gene identifiers")
  }
  if(!"gene" %in% names(stable.frame)){
    stop("Elastic-net stable-feature table does not contain a gene column")
  }
  selected.features = as.character(stable.frame\$gene)
  if(anyNA(selected.features) || anyDuplicated(selected.features)){
    stop("Elastic-net stable-feature table contains invalid gene identifiers")
  }
  missing.truth = setdiff(truedegs, allgenes)
  if(length(missing.truth) > 0){
    stop(sprintf("Synthetic truth set contains %d untested genes", length(missing.truth)))
  }
  missing.selected = setdiff(selected.features, allgenes)
  if(length(missing.selected) > 0){
    stop(sprintf("Elastic net selected %d genes outside the tested universe", length(missing.selected)))
  }
  n.tested = length(allgenes)
  n.truth = length(truedegs)
  n.background = n.tested - n.truth
  if(n.truth == 0) stop("Synthetic benchmark contains no true DEGs")
  if(n.background == 0) stop("Synthetic benchmark contains no background genes")

  # Classify every gene in the fixed synthetic truth universe
  truth = allgenes %in% truedegs
  called = allgenes %in% selected.features
  true.positives = sum(called & truth)
  false.positives = sum(called & !truth)
  true.negatives = sum(!called & !truth)
  false.negatives = sum(!called & truth)
  n.called = sum(called)

  # Calculate the common synthetic benchmark metrics
  power = true.positives / n.truth
  FDP = if(n.called == 0) 0 else false.positives / n.called
  youden.index = (true.positives / n.truth) -
                 (false.positives / n.background)
  metric.status = if(n.called == 0) "no_calls" else "ok"

  # Save the exact common full-synthetic result schema
  outframe = data.frame(
    method = "Stability-selected elastic-net features",
    power = power,
    FDP = FDP,
    youden_index = youden.index,
    true_positives = true.positives,
    false_positives = false.positives,
    true_negatives = true.negatives,
    false_negatives = false.negatives,
    n_tested = n.tested,
    n_called = n.called,
    n_truth = n.truth,
    n_background = n.background,
    metric_status = metric.status,
    replicate = as.integer($x),
    refnum = as.integer($refnum),
    altnum = as.integer($altnum),
    totalnum = as.integer($refnum) + as.integer($altnum),
    depth = as.numeric($depth)
  )
  write.csv(outframe, "outframe.csv", row.names = FALSE)
  RSCRIPT
  """
}
