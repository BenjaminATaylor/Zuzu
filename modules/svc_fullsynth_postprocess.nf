process SVC_FULLSYNTH_POSTPROCESS{

    debug false

    input:
    tuple val(depth), 
          path(truedegs), //"trueDEGs.RData"
          val(refnum),
          path(inframe),
          path(plotframe) // Included as input for ease of naming, but not used here 

    output:
    path("outframe.csv")

    script:
    """
    #!/usr/bin/env Rscript
    library("tidyverse", quietly = TRUE)
    library("edgeR", quietly = TRUE)
    library("pROC", quietly = TRUE)

    samplenum = $refnum
    truedegs = get(load("$truedegs"))
    depth = $depth
    inframe = read.csv("$inframe", row.names = 1)

    ## Generate metrics 

    # 1. Power: the proportion of treu DEGs identified as DEGs in this comparison
    degs = row.names(subset(inframe,DEGstatus==1))
    power = length(which(truedegs %in% degs))/length(truedegs)

    # 2. FDP: the proportion of all DEGs that are false positives
    FDP = 1-(length(which(degs %in% truedegs))/length(degs))

    # 3. ROC AUC
    allgenes = row.names(inframe)
    truelabels = as.numeric(allgenes %in% truedegs)
    labels = as.numeric(allgenes %in% degs)
    AUC = auc(truelabels, labels)
    # Remember that 0.5 = a truly random ROC, so for best effect we do (0.5-AUC)*2
    normAUC = (AUC-0.5)*2 # Greater deviation from 0 = better discrimination

    # Save outputs
    outframe = data.frame(method = "SVC", power = power, FDP = FDP, normAUC = normAUC, samplenum = samplenum, depth = depth)
    write.csv(outframe, file="outframe.csv", row.names = FALSE)

    """
}