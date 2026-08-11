process EDGER_QUASI {

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
    library("edgeR", quietly = TRUE)

    DEBUG=FALSE

    samplesheet = read.csv("$samplesheet")
    quasi.frame = read.csv("$quasiframe", row.names = 1, check.names = F)
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

    true.positives = NA_integer_
    false.positives = NA_integer_
    true.negatives = NA_integer_
    false.negatives = NA_integer_
    n.called = NA_integer_
    metric.status = "insufficient_truth"

    # Generate output metrics only if there are more than 10 'true' DEGs, otherwise output NAs
    if(length(truedegs)>10){

        #edgeR analysis using quasi-permuted data
        dds.edge = DGEList(counts=quasi.frame,group=samplesheet\$phenotype)
        dds.edge = normLibSizes(dds.edge)
        design.edge = model.matrix(~phenotype, data = samplesheet)
        phenotype.coef = which(attr(design.edge, "assign") == 1L)
        stopifnot(length(phenotype.coef) == 1L)
        stopifnot(all(design.edge[samplesheet\$phenotype == reflevel, phenotype.coef] == 0),
                  all(design.edge[samplesheet\$phenotype == altlevel, phenotype.coef] == 1))
        edge.contrast = rep(0, ncol(design.edge))
        edge.contrast[phenotype.coef] = 1
        dds.edge.est = estimateDisp(dds.edge,design.edge)
        edge.fit = glmQLFit(dds.edge.est,design.edge)
        edge.qlf = glmQLFTest(edge.fit,contrast=edge.contrast)
        edge.out.quasi = edge.qlf\$table %>% mutate(padj = p.adjust(PValue, method = "BH")) 
        stopifnot(identical(row.names(edge.out.quasi), allgenes))

        # Classify every gene in the fixed quasi-permutation universe.
        quasidegs = row.names(subset(edge.out.quasi,padj<0.05))
        truth = allgenes %in% truedegs
        called = allgenes %in% quasidegs
        true.positives = sum(called & truth)
        false.positives = sum(called & !truth)
        true.negatives = sum(!called & !truth)
        false.negatives = sum(!called & truth)
        n.called = sum(called)

        # 1. Power: the proportion of retained anchor genes identified as DEGs in this comparison
        power = true.positives / n.truth

        # 2. FDP: the proportion of all identified DEGs that fall outside the retained anchor set
        # FDP is defined as zero if no DEGs are identified
        FDP = if(n.called == 0) 0 else false.positives / n.called

        # 3. Youden's index: the anchor-recovery rate minus the background-call rate at the DEG threshold
        if(n.background == 0){
            youden_index = NA_real_
            metric.status = "youden_undefined_no_background"
        } else {
            youden_index = (true.positives / n.truth) -
                           (false.positives / n.background)
            metric.status = if(n.called == 0) "no_calls" else "ok"
        }

    } else {
        power = NA_real_
        FDP = NA_real_
        youden_index = NA_real_
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
