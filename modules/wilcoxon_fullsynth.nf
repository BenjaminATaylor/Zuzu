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

    # 3. Youden's index at the DEG threshold
    allgenes = row.names(wilcox.out)
    truth = allgenes %in% truedegs
    called = allgenes %in% degs
    sensitivity = sum(called & truth) / sum(truth)
    false_positive_rate = sum(called & !truth) / sum(!truth)
    youden_index = sensitivity - false_positive_rate

    # Save outputs
    outframe = data.frame(method = "Wilcoxon", power = power, FDP = FDP, youden_index = youden_index, samplenum = samplenum, depth = depth)
    write.csv(outframe, file="outframe.csv", row.names = FALSE)

    """
}
