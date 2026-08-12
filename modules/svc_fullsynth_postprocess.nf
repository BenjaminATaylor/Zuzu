process SVC_FULLSYNTH_POSTPROCESS{

    debug false

    input:
    tuple val(depth), 
          path(truedegs), //"trueDEGs.RData"
          val(refnum),
          val(altnum),
          val(x),
          path(inframe)

    output:
    path("outframe.csv")

    script:
    """
    #!/usr/bin/env Rscript
    library("tidyverse", quietly = TRUE)
    library("edgeR", quietly = TRUE)

    ref.num = as.integer($refnum)
    alt.num = as.integer($altnum)
    total.num = ref.num + alt.num
    replicate.id = $x
    truedegs = get(load("$truedegs"))
    depth = $depth
    inframe = read.csv("$inframe", row.names = 1)

    # Validate the requested group sizes
    stopifnot(ref.num > 0L, alt.num > 0L)

    # Validate the fixed synthetic evaluation universe
    allgenes = row.names(inframe)
    if(anyDuplicated(allgenes)) stop("SVC output contains duplicate gene identifiers")
    if(anyDuplicated(truedegs)) stop("Synthetic truth set contains duplicate gene identifiers")
    missing.truth = setdiff(truedegs, allgenes)
    if(length(missing.truth) > 0){
      stop(sprintf("Synthetic truth set contains %d untested genes", length(missing.truth)))
    }
    n.tested = length(allgenes)
    n.truth = length(truedegs)
    n.background = n.tested - n.truth
    if(n.truth == 0) stop("Synthetic benchmark contains no true DEGs")
    if(n.background == 0) stop("Synthetic benchmark contains no background genes")

    ## Generate metrics 

    # Classify every gene in the fixed synthetic truth universe.
    selected.features = row.names(subset(inframe, selected_feature == 1))
    truth = allgenes %in% truedegs
    called = allgenes %in% selected.features
    true.positives = sum(called & truth)
    false.positives = sum(called & !truth)
    true.negatives = sum(!called & !truth)
    false.negatives = sum(!called & truth)
    n.called = sum(called)

    # 1. Power: the proportion of true DEGs selected as features
    power = true.positives / n.truth

    # 2. FDP: the proportion of all selected features that are false positives
    # FDP is defined as zero if no features are selected
    FDP = if(n.called == 0) 0 else false.positives / n.called

    # 3. Youden's index: the true-positive rate minus the false-positive rate at the selected-feature threshold
    youden_index = (true.positives / n.truth) -
                   (false.positives / n.background)
    metric.status = if(n.called == 0) "no_calls" else "ok"

    # Save outputs
    outframe = data.frame(
      method = "SVC-RFECV candidates",
      power = power,
      FDP = FDP,
      youden_index = youden_index,
      true_positives = true.positives,
      false_positives = false.positives,
      true_negatives = true.negatives,
      false_negatives = false.negatives,
      n_tested = n.tested,
      n_called = n.called,
      n_truth = n.truth,
      n_background = n.background,
      metric_status = metric.status,
      replicate = replicate.id,
      refnum = ref.num,
      altnum = alt.num,
      totalnum = total.num,
      depth = depth
    )
    write.csv(outframe, file="outframe.csv", row.names = FALSE)

    """
}
