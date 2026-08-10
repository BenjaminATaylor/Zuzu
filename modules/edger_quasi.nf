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
    library("pROC", quietly = TRUE)

    DEBUG=FALSE

    samplesheet = read.csv("$samplesheet")
    quasi.frame = read.csv("$quasiframe", row.names = 1, check.names = F)
    truedegs = get(load("$truedegs"))
    replicate.id = $x
    ref.num = as.integer($refnum)
    alt.num = as.integer($altnum)
    reflevel = "$params.reflevel"

    # The quasi sample sheet has crossed another CSV boundary; restore its levels.
    samplesheet\$phenotype = relevel(factor(samplesheet\$phenotype), ref = reflevel)
    stopifnot(nlevels(samplesheet\$phenotype) == 2L,
              identical(levels(samplesheet\$phenotype)[1], reflevel))
    altlevel = setdiff(levels(samplesheet\$phenotype), reflevel)
    stopifnot(length(altlevel) == 1L)

    observed.refnum = sum(samplesheet\$phenotype == reflevel)
    observed.altnum = sum(samplesheet\$phenotype != reflevel)
    stopifnot(observed.refnum == ref.num)
    stopifnot(observed.altnum == alt.num)
    stopifnot(nrow(samplesheet) == ref.num + alt.num)
    stopifnot(ncol(quasi.frame) == ref.num + alt.num)
    stopifnot(identical(colnames(quasi.frame), samplesheet\$sample))

    # Generate output metrics only if there are at least 10 'true' DEGs, otherwise outputs NAs
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

    } else {
        power = NA
        FDP = NA
        normAUC = NA
    }

    outframe = data.frame(power = power,
                          FDP = FDP,
                          normAUC = normAUC,
                          replicate = replicate.id,
                          refnum = ref.num,
                          altnum = alt.num,
                          totalnum = nrow(samplesheet))
    write.csv(outframe, file="outframe.csv", row.names = FALSE)
"""

}
