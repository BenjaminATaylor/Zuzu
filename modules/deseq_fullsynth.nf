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
    library("pROC", quietly = TRUE)

    samplesheet = read.csv("$samplesheet")
    countsframe = read.csv("$countsframe", row.names = 1, check.names = F)
    samplenum = $refnum
    truedegs = get(load("$truedegs"))
    depth = $depth

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

    ## Generate metrics 

    # 1. Power: the proportion of treu DEGs identified as DEGs in this comparison
    degs = row.names(subset(deseq.results,padj<0.05))
    power = length(which(truedegs %in% degs))/length(truedegs)

    # 2. FDP: the proportion of all DEGs that are false positives
    FDP = 1-(length(which(degs %in% truedegs))/length(degs))

    # 3. ROC AUC
    allgenes = row.names(deseq.results)
    truelabels = as.numeric(allgenes %in% truedegs)
    labels = as.numeric(allgenes %in% degs)
    AUC = auc(truelabels, labels)
    # Remember that 0.5 = a truly random ROC, so for best effect we do (0.5-AUC)*2
    normAUC = (AUC-0.5)*2 # Greater deviation from 0 = better discrimination

    # Save outputs
    outframe = data.frame(method = "DESeq2", power = power, FDP = FDP, normAUC = normAUC, samplenum = samplenum, depth = depth)
    write.csv(outframe, file="outframe.csv", row.names = FALSE)

    """
}
