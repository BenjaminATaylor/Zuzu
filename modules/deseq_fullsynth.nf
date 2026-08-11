process DESEQ_FULLSYNTH{

    label 'mid_job'

    input:
    tuple val(depth), 
          path(truedegs), //"trueDEGs.RData"
          path(samplesheet), //"synthsheet.csv"
          path(countsframe), //"synthcounts.csv"
          val(refnum),
          val(x)

    output:
    path("outframe.csv")

    script:
    """
    #!/usr/bin/env Rscript
    library("tidyverse", quietly = TRUE)
    library("DESeq2", quietly = TRUE)

    samplesheet = read.csv("$samplesheet")
    countsframe = read.csv("$countsframe", row.names = 1, check.names = F)
    samplenum = $refnum
    replicate.id = $x
    truedegs = get(load("$truedegs"))
    depth = $depth

    # Validate the fixed synthetic evaluation universe
    allgenes = row.names(countsframe)
    if(anyDuplicated(allgenes)) stop("Synthetic count matrix contains duplicate gene identifiers")
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

    # Restore sample sheet factor levels
    reflevel = "$params.reflevel"
    samplesheet\$phenotype = relevel(factor(samplesheet\$phenotype), ref = reflevel)
    stopifnot(nlevels(samplesheet\$phenotype) == 2L,
              identical(levels(samplesheet\$phenotype)[1], reflevel))
    altlevel = setdiff(levels(samplesheet\$phenotype), reflevel)
    stopifnot(length(altlevel) == 1L)

    # Run DESeq for this dataset
    dds = DESeqDataSetFromMatrix(countData = countsframe,
                                colData = samplesheet,
                                design = as.formula(~phenotype))
    # Run the default analysis for DESeq2
    dds.deg = DESeq(dds, fitType = "parametric", betaPrior = FALSE)
    deseq.results = results(dds.deg,
                            contrast = c("phenotype", altlevel, reflevel))
    stopifnot(identical(row.names(deseq.results), allgenes))

    ## Generate metrics 

    # Classify every gene in the fixed synthetic truth universe.
    degs = row.names(subset(deseq.results,padj<0.05))
    truth = allgenes %in% truedegs
    called = allgenes %in% degs
    true.positives = sum(called & truth)
    false.positives = sum(called & !truth)
    true.negatives = sum(!called & !truth)
    false.negatives = sum(!called & truth)
    n.called = sum(called)

    # 1. Power: the proportion of true DEGs identified as DEGs in this comparison
    power = true.positives / n.truth

    # 2. FDP: the proportion of all identified DEGs that are false positives
    # FDP is defined as zero if no DEGs are identified
    FDP = if(n.called == 0) 0 else false.positives / n.called

    # 3. Youden's index: the true-positive rate minus the false-positive rate at the DEG threshold
    youden_index = (true.positives / n.truth) -
                   (false.positives / n.background)
    metric.status = if(n.called == 0) "no_calls" else "ok"

    # Save outputs
    outframe = data.frame(
      method = "DESeq2",
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
      samplenum = samplenum,
      depth = depth
    )
    write.csv(outframe, file="outframe.csv", row.names = FALSE)

    """
}
