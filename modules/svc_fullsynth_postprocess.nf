process SVC_FULLSYNTH_POSTPROCESS{

    debug false

    input:
    tuple val(depth), 
          path(truedegs), //"trueDEGs.RData"
          val(refnum),
          path(inframe),
          path(plotframe), // Included as input for ease of naming, but not used here
          path(pickle) // Included as input for ease of naming, but not used here

    output:
    path("outframe.csv")

    script:
    """
    #!/usr/bin/env Rscript
    library("tidyverse", quietly = TRUE)
    library("edgeR", quietly = TRUE)

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

    # 3. Youden's index at the selected-feature threshold
    allgenes = row.names(inframe)
    truth = allgenes %in% truedegs
    called = allgenes %in% degs
    sensitivity = sum(called & truth) / sum(truth)
    false_positive_rate = sum(called & !truth) / sum(!truth)
    youden_index = sensitivity - false_positive_rate

    # Save outputs
    outframe = data.frame(method = "SVC", power = power, FDP = FDP, youden_index = youden_index, samplenum = samplenum, depth = depth)
    write.csv(outframe, file="outframe.csv", row.names = FALSE)

    """
}
