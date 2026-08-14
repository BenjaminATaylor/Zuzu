process INTER_OBSERVER_BASIC{

  debug false

  publishDir "$params.outdir"
  label 'mid_job'

  input:
  path deseq_table
  path deseq_dds
  path edger_table
  path wilcox_table
  path SVC_table, stageAs: 'svc_basic.csv'
  path enet_table, stageAs: 'enet_basic.csv'

  output:
  path "upSet_basic.pdf"
  path "basic_overlap_violins.pdf"
  path "deseq_poorfits.RData", emit: deseq_poorfits

    script:
    """
    #!/usr/bin/env Rscript
    library("tidyverse", quietly = TRUE)
    library("UpSetR", quietly = TRUE)
    library("DESeq2")
    library("vcd")
    library("ggpubr")

    deseq.full = get(load("$deseq_dds"))

    #hypergeometric overlap test function
    gene_overlap_test = function(list1, list2, background,verbose=TRUE,lower.tail=FALSE){
        n_A = length(list1)
        n_B = length(list2)
        n_C = length(background)
        n_A_B = length(intersect(list1,list2))
        hyp = phyper(n_A_B - 1, n_A, n_C-n_A, n_B, lower.tail = lower.tail)
        expc = round((n_A/n_C)*n_B)
        jac = n_A_B/(n_A+n_B-n_A_B)
        if(verbose){print(paste0(n_A_B," matching from lists of lengths ",n_A," and ",
                                n_B,"; p=",round(hyp,4),". Expected intersection length = ",expc))}
        res = list(hypergeom = round(hyp,4),
                    jaccard = jac,
                    intersect = n_A_B)
        return(res)
    }
    ## DESeq inputs
    # DEGs from full model
    deseq.table = read.csv("$deseq_table", row.names = 1)
    edger.table = read.csv("$edger_table", row.names = 1)
    wilcox.table = read.csv("$wilcox_table", row.names = 1)
    getdegs = function(df, ml=FALSE){ 
    if(ml){
            return(row.names(subset(df, selected_feature == 1)))
        } else { 
            return(row.names(subset(df, padj<0.05))) 
        } 
    }
    degs.deseq = getdegs(deseq.table)
    degs.edger = getdegs(edger.table)
    degs.wilcox = getdegs(wilcox.table)
    # Only include ML outputs if specified by user
    # TODO: Fix this to work for ML outputs also
    if("false" == "true"){
        svc.table = read.csv("$SVC_table", row.names = 1)
        degs.svc = getdegs(svc.table, ml=TRUE)
    }

    # Combine adjusted p-values across the shared gene order
    compframe = cbind(select(deseq.table,padj), 
                    select(edger.table,padj),
                    select(wilcox.table,padj))  %>% 
    remove_rownames() %>%
    `colnames<-`(c("deseq","edger","wilcox"))

    pairs = t(combn(colnames(compframe),2))
            
    for(i in 1:nrow(pairs)){
  
        a = compframe[,pairs[i,1]]
        b = compframe[,pairs[i,2]]
        
        # First let's test the correlations of the p-values of the two sets (only relevant where we have pvals)
        # RRHO would be nice here but neither of the two existing packages (RRHO and RedRibbon) are suitable- the former due to bugs and the latter due to instability
        # TODO: apply my own RRHO algorithm
        cor.kendall = cor.test(a, b, 
                                alternative = "two.sided",
                                method = "kendall")\$estimate
        cor.pearson = cor.test(a, b, 
                                alternative = "two.sided",
                                method = "pearson")\$estimate
        plot(a,b)
        
        # As expected, the pval correlations themselves are usually very high. The more interesting question is how similar each pair of methods is in terms of the specific set of DEGs identified
        
        alldegs = row.names(compframe)
        a = get(paste0("degs.",pairs[i,1]))
        b = get(paste0("degs.",pairs[i,2]))
        
        thisgot = gene_overlap_test(a, b, alldegs, verbose = FALSE)
        jacc = thisgot\$jaccard
        hyprp = thisgot\$hypergeom
    }

    # Also generate an UpSet showing the baseline gene-list intersections
    listInput.combined <- list("DESeq2" = degs.deseq,
                            "edgeR" = degs.edger,
                            "Wilcoxon" = degs.wilcox)

    upset.lists = listInput.combined
    if("$params.enetstep" == "true"){
        enet.table = read.csv("$enet_table", row.names = 1, check.names = FALSE)
        degs.enet = row.names(subset(enet.table, stable_feature))
        upset.lists[["Elastic net"]] = degs.enet
    }

    print(upset.combined <- upset(fromList(upset.lists), order.by = "degree",
                                point.size = 4, line.size = 2,
                                mainbar.y.label = "Gene-list intersections",
                                sets.x.label = "Selected features",
                                empty.intersections = "on", text.scale = 3.2)
    )

    pdf(file="upSet_basic.pdf",width = 24, height = 18) ; upset.combined ;dev.off()    
        
    # From Li et al.: we can estimate each gene's poorness-of-fit to the nbinom distribution under each model
    # Doing so let's us asked whether mismatched genes are disproportionately likely to fail this assumption
    # NB to fit li et al we should use the package-normalized methods for each, but for testing we start with just the DESeq2 fits
    deseq.dispersions = dispersions(deseq.full)
    deseq.normcounts = as.list(data.frame(t(counts(deseq.full, normalized=TRUE))))

    getfit = function(counts,disp){
    
    thesecounts = round(counts)
    fit = summary(goodfit(thesecounts, type = "nbinomial", 
                    par = list(size = 1/disp)))
    return(-log10(fit[3]))
    
    }

    deseq.poorfits = mapply(getfit, deseq.normcounts, deseq.dispersions)
    save(deseq.poorfits,file = "deseq_poorfits.RData") # Save for downstream

    # Divide genes into sets based on how many baseline gene lists selected them
    # Get IDs used by the expanded UpSet
    thesegenes = unlist(upset.lists, use.names = FALSE)
    thesegenes = thesegenes[!duplicated(thesegenes)]
    # Count membership across the conventional and elastic-net lists
    gene.list.membership = vapply(
      upset.lists,
      function(degs) thesegenes %in% degs,
      logical(length(thesegenes))
    )
    intersects = data.frame(gene = thesegenes, sum = rowSums(gene.list.membership))
    ggframe = left_join(intersects, rownames_to_column(data.frame(deseq.poorfits),var = "gene")) %>%
    `colnames<-`(c("gene","sum","poorfit"))

    # Compare the negative-binomial fit between every observed membership count
    observed.counts = sort(unique(as.character(ggframe\$sum)))
    poorfit.comparisons = if(length(observed.counts) > 1){
      combn(observed.counts, 2, simplify = FALSE)
    } else {
      list()
    }

    gg.poor = ggplot(ggframe, aes(x = as.factor(sum), y = poorfit)) +
    geom_violin() + 
    labs(x = "Number of baseline gene lists selecting each gene",
        y = "Poorness of fit to Negative Binomial distribution") +
    theme_bw()

    if(length(poorfit.comparisons) > 0){
      gg.poor = gg.poor +
        stat_compare_means(method = "wilcox.test",
                           label = "p.signif",
                           size = 5,
                           comparisons = poorfit.comparisons)
    }

    ggsave(gg.poor, 
        filename = "basic_overlap_violins.pdf",
        device = "pdf", bg = "transparent",
        width =  12, height = 12, units = "cm")
    """

}
