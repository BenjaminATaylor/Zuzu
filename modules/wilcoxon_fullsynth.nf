process WILCOXON_FULLSYNTH{

    input:
    tuple val(depth), 
          path(truedegs), //"trueDEGs.RData"
          path(samplesheet), //"synthsheet.csv"
          path(countsframe), //"synthcounts.csv"
          val(refnum),
          val(altnum),
          val(x)

    output:
    path("outframe.csv")

    script:
    """
    #!/usr/bin/env Rscript
    library("tidyverse", quietly = TRUE)
    library("edgeR", quietly = TRUE)

    samplesheet = read.csv("$samplesheet")
    countsframe = read.csv("$countsframe", row.names = 1, check.names = F)
    ref.num = as.integer($refnum)
    alt.num = as.integer($altnum)
    total.num = ref.num + alt.num
    replicate.id = $x
    truedegs = get(load("$truedegs"))
    depth = $depth
    reflevel = "$params.reflevel"

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

    # Validate the requested group sizes and sample alignment
    stopifnot(ref.num > 0L,
              alt.num > 0L,
              nrow(samplesheet) == total.num,
              ncol(countsframe) == total.num,
              !anyNA(samplesheet\$sample),
              !anyNA(samplesheet\$phenotype),
              !anyDuplicated(samplesheet\$sample),
              !anyDuplicated(colnames(countsframe)),
              identical(colnames(countsframe), as.character(samplesheet\$sample)),
              length(unique(samplesheet\$phenotype)) == 2L,
              reflevel %in% samplesheet\$phenotype,
              sum(samplesheet\$phenotype == reflevel) == ref.num,
              sum(samplesheet\$phenotype != reflevel) == alt.num)

    #Wilcox analysis for this dataset
    countsframe.dge = DGEList(counts=countsframe,group=samplesheet\$phenotype)
    countsframe.dge = calcNormFactors(countsframe.dge,method="TMM")
    countsframe.dge = cpm(countsframe.dge) %>% as.data.frame()
    stopifnot(identical(row.names(countsframe.dge), allgenes))

    # Split data into ref and alt counts
    refcounts = select(countsframe.dge, subset(samplesheet, phenotype == reflevel)\$sample)
    altcounts = select(countsframe.dge, subset(samplesheet, phenotype != reflevel)\$sample)
    # Run wilcox tests on each gene
    gene.wilcox = function(x){
        this.wilcox = wilcox.test(as.numeric(refcounts[x,]), as.numeric(altcounts[x,]))
        return(this.wilcox\$p.value)
    }
    wilcox.p = sapply(allgenes, gene.wilcox)
    # Adjust to FDR
    wilcox.out = p.adjust(wilcox.p,method = "BH") %>% 
        `names<-`(allgenes) %>% 
        data.frame() %>%
        `colnames<-`("padj")

    ## Generate metrics 

    # Classify every gene in the fixed synthetic truth universe.
    degs = row.names(subset(wilcox.out,padj<0.05))
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
      method = "Wilcoxon",
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
