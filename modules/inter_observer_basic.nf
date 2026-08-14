process INTER_OBSERVER_BASIC{

  debug false

  publishDir "$params.outdir", pattern: '*.pdf'
  publishDir "$params.outdir", pattern: '*.RData'
  publishDir "$params.outdir", pattern: 'baseline_*.csv', mode: 'copy'
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
  path "baseline_membership.csv", emit: membership
  path "baseline_pairwise_overlap.csv"
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

    ## DESeq inputs
    # DEGs from full model
    deseq.table = read.csv("$deseq_table", row.names = 1, check.names = FALSE)
    edger.table = read.csv("$edger_table", row.names = 1, check.names = FALSE)
    wilcox.table = read.csv("$wilcox_table", row.names = 1, check.names = FALSE)
    analysis.genes = row.names(deseq.table)
    if(!setequal(analysis.genes, row.names(edger.table)) ||
       !setequal(analysis.genes, row.names(wilcox.table))){
        stop("Conventional baseline tables do not share the same gene universe")
    }
    edger.table = edger.table[analysis.genes,,drop = FALSE]
    wilcox.table = wilcox.table[analysis.genes,,drop = FALSE]

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

    # Also generate an UpSet showing the baseline gene-list intersections
    listInput.combined <- list("DESeq2" = degs.deseq,
                            "edgeR" = degs.edger,
                            "Wilcoxon" = degs.wilcox)

    upset.lists = listInput.combined
    if("$params.enetstep" == "true"){
        enet.table = read.csv("$enet_table", row.names = 1, check.names = FALSE)
        if(!setequal(analysis.genes, row.names(enet.table))){
            stop("Elastic-net and conventional baselines do not share the same gene universe")
        }
        enet.table = enet.table[analysis.genes,,drop = FALSE]
        degs.enet = row.names(subset(enet.table, stable_feature))
        upset.lists[["Elastic net"]] = degs.enet
    }

    # Keep adjusted p-value correlations for pairs where both methods report them
    # RRHO would be nice here but neither existing package was sufficiently stable
    # TODO: apply my own RRHO algorithm
    conventional.padj = list(
      "DESeq2" = deseq.table\$padj,
      "edgeR" = edger.table\$padj,
      "Wilcoxon" = wilcox.table\$padj
    )
    get.padj.correlation = function(method.pair, correlation.method){
      if(!all(method.pair %in% names(conventional.padj))){
        return(NA_real_)
      }
      a = conventional.padj[[method.pair[1]]]
      b = conventional.padj[[method.pair[2]]]
      complete = complete.cases(a, b)
      if(sum(complete) < 2L || sd(a[complete]) == 0 || sd(b[complete]) == 0){
        return(NA_real_)
      }
      cor(a[complete], b[complete], method = correlation.method)
    }

    # Retain every gene's membership in each enabled baseline list
    membership.matrix = vapply(
      upset.lists,
      function(degs) analysis.genes %in% degs,
      logical(length(analysis.genes))
    )
    baseline.membership = data.frame(
      gene = analysis.genes,
      membership.matrix,
      check.names = FALSE
    )
    if(!"Elastic net" %in% names(baseline.membership)){
      baseline.membership[["Elastic net"]] = NA
    }
    baseline.membership[["n_methods_run"]] = length(upset.lists)
    baseline.membership[["n_methods_selected"]] = rowSums(membership.matrix)
    write.csv(baseline.membership, "baseline_membership.csv", row.names = FALSE)

    # Retain pairwise gene-list agreement over the shared tested universe
    method.pairs = combn(names(upset.lists), 2, simplify = FALSE)
    baseline.pairwise.overlap = lapply(method.pairs, function(method.pair){
      list.1 = unique(upset.lists[[method.pair[1]]])
      list.2 = unique(upset.lists[[method.pair[2]]])
      intersection.size = length(intersect(list.1, list.2))
      union.size = length(union(list.1, list.2))
      data.frame(
        method_1 = method.pair[1],
        method_2 = method.pair[2],
        n_tested = length(analysis.genes),
        n_selected_1 = length(list.1),
        n_selected_2 = length(list.2),
        selection_proportion_1 = length(list.1) / length(analysis.genes),
        selection_proportion_2 = length(list.2) / length(analysis.genes),
        intersection = intersection.size,
        union = union.size,
        jaccard = if(union.size == 0L) NA_real_ else intersection.size / union.size,
        expected_intersection = length(list.1) * length(list.2) / length(analysis.genes),
        hypergeometric_p_value = phyper(
          intersection.size - 1,
          length(list.1),
          length(analysis.genes) - length(list.1),
          length(list.2),
          lower.tail = FALSE
        ),
        padj_kendall_correlation = get.padj.correlation(method.pair, "kendall"),
        padj_pearson_correlation = get.padj.correlation(method.pair, "pearson"),
        overlap_status = if(union.size == 0L) "both_lists_empty" else "ok",
        stringsAsFactors = FALSE
      )
    }) %>% bind_rows()
    write.csv(
      baseline.pairwise.overlap,
      "baseline_pairwise_overlap.csv",
      row.names = FALSE
    )

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
