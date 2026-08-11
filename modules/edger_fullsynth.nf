process EDGER_FULLSYNTH{

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
    library("edgeR", quietly = TRUE)

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

    # Run edgeR for this dataset
    dds.edge = DGEList(counts=countsframe,group=samplesheet\$phenotype)
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
    edge.out = edge.qlf\$table %>% mutate(padj = p.adjust(PValue, method = "BH")) 

    ## Generate metrics 

    # 1. Power: the proportion of treu DEGs identified as DEGs in this comparison
    degs = row.names(subset(edge.out,padj<0.05))
    power = length(which(truedegs %in% degs))/length(truedegs)

    # 2. FDP: the proportion of all DEGs that are false positives
    FDP = 1-(length(which(degs %in% truedegs))/length(degs))

    # 3. Youden's index at the DEG threshold
    allgenes = row.names(edge.out)
    truth = allgenes %in% truedegs
    called = allgenes %in% degs
    sensitivity = sum(called & truth) / sum(truth)
    false_positive_rate = sum(called & !truth) / sum(!truth)
    youden_index = sensitivity - false_positive_rate

    # Save outputs
    outframe = data.frame(method = "edgeR", power = power, FDP = FDP, youden_index = youden_index, samplenum = samplenum, depth = depth)
    write.csv(outframe, file="outframe.csv", row.names = FALSE)

    """
}
