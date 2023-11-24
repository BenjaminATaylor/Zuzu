process EDGER_FULLSYNTH{

    input:
    tuple path(depth), //"depth.RData"
          path(truedegs), //"trueDEGs.RData"
          path(samplesheet), //"synthsheet.csv"
          path(countsframe), //"synthcounts.csv"
          val(refnum)

    output:
    path("outframe.csv")

    script:
    """
    #!/usr/bin/env Rscript
    library("tidyverse", quietly = TRUE)
    library("edgeR", quietly = TRUE)
    library("pROC", quietly = TRUE)

    samplesheet = read.csv("$samplesheet")
    countsframe = read.csv("$countsframe", row.names = 1, check.names = F)
    samplenum = $refnum
    truedegs = get(load("$truedegs"))
    depth = get(load("$depth"))

    # Run edgeR for this dataset
    dds.edge = DGEList(counts=countsframe,group=samplesheet\$phenotype)
    dds.edge = normLibSizes(dds.edge)
    design.edge = model.matrix(~samplesheet\$phenotype)
    dds.edge.est = estimateDisp(dds.edge,design.edge)
    edge.fit = glmQLFit(dds.edge.est,design.edge)
    edge.qlf = glmQLFTest(edge.fit,coef=2)
    edge.out = edge.qlf\$table %>% mutate(padj = p.adjust(PValue, method = "BH")) 

    ## Generate metrics 

    # 1. Power: the proportion of treu DEGs identified as DEGs in this comparison
    degs = row.names(subset(edge.out,padj<0.05))
    power = length(which(truedegs %in% degs))/length(truedegs)

    # 2. FDP: the proportion of all DEGs that are false positives
    FDP = 1-(length(which(degs %in% truedegs))/length(degs))

    # 3. ROC AUC
    allgenes = row.names(edge.out)
    truelabels = as.numeric(allgenes %in% truedegs)
    labels = as.numeric(allgenes %in% degs)
    AUC = auc(truelabels, labels)
    # Remember that 0.5 = a truly random ROC, so for best effect we do (0.5-AUC)*2
    normAUC = (AUC-0.5)*2 # Greater deviation from 0 = better discrimination

    # Save outputs
    outframe = data.frame(method = "edgeR", power = power, FDP = FDP, normAUC = normAUC, samplenum = samplenum, depth = depth)
    write.csv(outframe, file="outframe.csv", row.names = FALSE)

    """
}