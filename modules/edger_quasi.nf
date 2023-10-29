process EDGER_QUASI {
  
    input: 
    tuple path(quasiframe), path(truedegs), path(samplesheet), val(samplenum)

    output:
    path 'outframe.csv'

    script:
    """
    #!/usr/bin/env Rscript
    library("tidyverse", quietly = TRUE)
    library("edgeR", quietly = TRUE)
    library("pROC", quietly = TRUE)

    DEBUG=FALSE

    samplesheet = read.csv("$samplesheet")
    quasi.frame = read.csv("$quasiframe", row.names = 1, check.names = F)
    truedegs = get(load("$truedegs"))
    samplenum = $samplenum

    # Generate output metrics only if there are at least 50 'true' DEGs, otherwise outputs NAs
    if(length(truedegs)>50){

        #edgeR analysis using quasi-permuted data
        dds.edge = DGEList(counts=quasi.frame,group=samplesheet\$phenotype)
        dds.edge = normLibSizes(dds.edge)
        design.edge = model.matrix(~samplesheet\$phenotype)
        dds.edge.est = estimateDisp(dds.edge,design.edge)
        edge.fit = glmQLFit(dds.edge.est,design.edge)
        edge.qlf = glmQLFTest(edge.fit,coef=2)
        edge.out.quasi = edge.qlf\$table %>% mutate(padj = p.adjust(PValue, method = "BH")) 

        # 1. Power: the proportion of treu DEGs identified as DEGs in this comparison
        quasidegs = row.names(subset(edge.out.quasi,padj<0.05))
        power = length(which(truedegs %in% quasidegs))/length(truedegs)

        # 2. FDP: the proportion of all DEGs that are false positives
        FDP = 1-(length(which(quasidegs %in% truedegs))/length(quasidegs))

        # 3. ROC AUC
        allgenes = row.names(quasi.frame)
        truelabels = as.numeric(allgenes %in% truedegs)
        quasilabels = as.numeric(allgenes %in% quasidegs)
        AUC = auc(truelabels, quasilabels)
        # Remember that 0.5 = a truly random ROC, so for best effect we do (0.5-AUC)*2
        normAUC = (AUC-0.5)*2 # Greater deviation from 0 = better discrimination

        # Save outputs
        outframe = data.frame(power = power, FDP = FDP, normAUC = normAUC, samplenum = samplenum)
        write.csv(outframe, file="outframe.csv", row.names = FALSE)
    } else {
        outframe = data.frame(power = NA, FDP = NA, normAUC = NA, samplenum = samplenum)
        write.csv(outframe, file="outframe.csv", row.names = FALSE)
    }

"""

}