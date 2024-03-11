process WILCOXON_FULLSYNTH{

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
    library("pROC", quietly = TRUE)

    samplesheet = read.csv("$samplesheet")
    countsframe = read.csv("$countsframe", row.names = 1, check.names = F)
    samplenum = $refnum
    truedegs = get(load("$truedegs"))
    depth = $depth
    reflevel = "$params.reflevel"

    #Wilcox analysis for this dataset
    countsframe.dge = DGEList(counts=countsframe,group=samplesheet\$phenotype)
    countsframe.dge = calcNormFactors(countsframe.dge,method="TMM")
    countsframe.dge = cpm(countsframe.dge) %>% as.data.frame()
    allgenes = row.names(countsframe.dge)

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

    # 1. Power: the proportion of treu DEGs identified as DEGs in this comparison
    degs = row.names(subset(wilcox.out,padj<0.05))
    power = length(which(truedegs %in% degs))/length(truedegs)

    # 2. FDP: the proportion of all DEGs that are false positives
    FDP = 1-(length(which(degs %in% truedegs))/length(degs))

    # 3. ROC AUC
    allgenes = row.names(wilcox.out)
    truelabels = as.numeric(allgenes %in% truedegs)
    labels = as.numeric(allgenes %in% degs)
    AUC = auc(truelabels, labels)
    # Remember that 0.5 = a truly random ROC, so for best effect we do (0.5-AUC)*2
    normAUC = (AUC-0.5)*2 # Greater deviation from 0 = better discrimination

    # Save outputs
    outframe = data.frame(method = "Wilcoxon", power = power, FDP = FDP, normAUC = normAUC, samplenum = samplenum, depth = depth)
    write.csv(outframe, file="outframe.csv", row.names = FALSE)

    """
}