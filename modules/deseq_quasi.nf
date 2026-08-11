process DESEQ_QUASI {
  
  label 'mid_job'
  tag "replicate=${x},ref=${refnum},alt=${altnum}"
  
  input: 
  tuple path(quasiframe),
        path(truedegs),
        path(samplesheet),
        val(x),
        val(refnum),
        val(altnum)

  output:
  path 'outframe.csv'

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse", quietly = TRUE)
  library("DESeq2", quietly = TRUE)

  DEBUG=FALSE

  samplesheet = read.csv("$samplesheet")
  quasi.frame = read.csv("$quasiframe", row.names = 1, check.names = FALSE)
  truedegs = get(load("$truedegs"))
  replicate.id = $x
  ref.num = as.integer($refnum)
  alt.num = as.integer($altnum)
  reflevel = "$params.reflevel"

  # Restore sample sheet factor levels
  samplesheet\$phenotype = relevel(factor(samplesheet\$phenotype), ref = reflevel)
  stopifnot(nlevels(samplesheet\$phenotype) == 2L,
            identical(levels(samplesheet\$phenotype)[1], reflevel))
  altlevel = setdiff(levels(samplesheet\$phenotype), reflevel)
  stopifnot(length(altlevel) == 1L)

  # Validate sample sizes and count-matrix alignment
  observed.refnum = sum(samplesheet\$phenotype == reflevel)
  observed.altnum = sum(samplesheet\$phenotype != reflevel)
  stopifnot(observed.refnum == ref.num)
  stopifnot(observed.altnum == alt.num)
  stopifnot(nrow(samplesheet) == ref.num + alt.num)
  stopifnot(ncol(quasi.frame) == ref.num + alt.num)
  stopifnot(identical(colnames(quasi.frame), samplesheet\$sample))

  # Validate the anchors against the fixed quasi-permutation gene universe
  allgenes = row.names(quasi.frame)
  if(anyDuplicated(allgenes)) stop("Quasi count matrix contains duplicate gene identifiers")
  if(anyDuplicated(truedegs)) stop("Quasi truth set contains duplicate gene identifiers")
  missing.truth = setdiff(truedegs, allgenes)
  if(length(missing.truth) > 0){
    stop(sprintf("Quasi truth set contains %d untested genes", length(missing.truth)))
  }
  n.tested = length(allgenes)
  n.truth = length(truedegs)
  n.background = n.tested - n.truth

  # Run the default analysis for DESeq2
  dds = DESeqDataSetFromMatrix(countData = quasi.frame,
                              colData = samplesheet,
                              design = as.formula(~phenotype))
  dds.deg.quasi = DESeq(dds, fitType = "parametric", betaPrior = FALSE)
  deseq.results = results(dds.deg.quasi,
                          contrast = c("phenotype", altlevel, reflevel))
  stopifnot(identical(row.names(deseq.results), allgenes))

  # Classify every gene in the fixed quasi-permutation universe
  quasidegs = row.names(subset(deseq.results,padj<0.05))
  if(DEBUG){quasidegs = sample(row.names(deseq.results),200)}
  #if(DEBUG){quasidegs = sample(truedegs, 70)}
  truth = allgenes %in% truedegs
  called = allgenes %in% quasidegs
  true.positives = sum(called & truth)
  false.positives = sum(called & !truth)
  true.negatives = sum(!called & !truth)
  false.negatives = sum(!called & truth)
  n.called = sum(called)

  # 1. Power: the proportion of retained anchor genes identified as DEGs in this comparison
  # Power is undefined if no anchors were retained
  power = if(n.truth == 0) NA_real_ else true.positives / n.truth

  # 2. FDP: the proportion of all identified DEGs that fall outside the retained anchor set
  # FDP is defined as zero if no DEGs are identified
  FDP = if(n.called == 0) 0 else false.positives / n.called

  # 3. Youden's index: the anchor-recovery rate minus the background-call rate at the DEG threshold
  # Youden's index requires both anchor and background genes
  youden_index = if(n.truth == 0 || n.background == 0){
    NA_real_
  } else {
    (true.positives / n.truth) - (false.positives / n.background)
  }

  # Record the primary condition affecting metric interpretation
  metric.status = if(n.truth == 0){
    "no_anchors"
  } else if(n.background == 0){
    "youden_undefined_no_background"
  } else if(n.truth <= 10){
    "low_anchor_count"
  } else if(n.called == 0){
    "no_calls"
  } else {
    "ok"
  }

  outframe = data.frame(power = power,
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
                        totalnum = nrow(samplesheet))
  write.csv(outframe, file="outframe.csv", row.names = FALSE)
  """

}
