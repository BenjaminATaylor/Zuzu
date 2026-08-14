process PERMUTE_HISTS{

  publishDir "$params.outdir"

  input:
  path(deseq_infiles,  stageAs: 'deseq/deseq_table_?.csv')
  path(edger_infiles,  stageAs: 'edger/edger_table_?.csv')
  path(wilcox_infiles, stageAs: 'wilcox/wilcox_table_?.csv')
  path(SVC_infiles,    stageAs: 'svc/svc_table_?.csv')
  path(enet_infiles,   stageAs: 'enet/enet_stable_features_?.csv')
  path deseq_poorfits

  output:
  path 'permute_hist.pdf'
  path 'permute_violins.pdf', optional: true

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse", quietly = TRUE)
  library("reshape2", quietly = TRUE)
  library("ggpubr")

  nperms = as.numeric("$params.nperms")

  # Count how often each gene is selected across permutations
  read_permdegs_r = function(dir){
    lapply(list.files(dir, full.names = TRUE), function(x){
      res <- read.csv(x, row.names = 1)
      data.frame(result = (res\$padj < 0.05), row.names = row.names(res))
    }) %>% list_cbind() %>% rowSums(na.rm = TRUE)
  }

  read_permfeatures_py = function(dir){
    lapply(list.files(dir, full.names = TRUE), function(x){
      res <- read.csv(x, row.names = 1)
      data.frame(result = (res\$selected_feature == 1), row.names = row.names(res))
    }) %>% list_cbind() %>% rowSums(na.rm = TRUE)
  }

  read_permfeatures_enet = function(dir, universe){
    result = setNames(integer(length(universe)), universe)
    files = list.files(dir, full.names = TRUE)
    if(length(files) == 0L){
      stop("No elastic-net permutation feature tables were supplied")
    }
    for(file in files){
      frame = read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
      if(!"gene" %in% names(frame)){
        stop(sprintf("Elastic-net permutation table lacks gene: %s", basename(file)))
      }
      genes = as.character(frame\$gene)
      if(anyNA(genes) || anyDuplicated(genes)){
        stop(sprintf("Elastic-net permutation table has invalid genes: %s", basename(file)))
      }
      missing = setdiff(genes, universe)
      if(length(missing) > 0L){
        stop(sprintf(
          "Elastic-net permutation table contains %d genes outside the common universe",
          length(missing)
        ))
      }
      result[genes] = result[genes] + 1L
    }
    result
  }

  deseq.permdegs  = read_permdegs_r("deseq")
  edger.permdegs  = read_permdegs_r("edger")
  wilcox.permdegs = read_permdegs_r("wilcox")

  # Combine all requested methods over the common gene order
  permhist.wide = cbind(deseq.permdegs,
                        edger.permdegs,
                        wilcox.permdegs) %>%
    `colnames<-`(c("DESeq2", "edgeR", "Wilcoxon"))

  if("$params.mlstep" == "true"){
    SVC.permfeatures = read_permfeatures_py("svc")
    if(!setequal(names(SVC.permfeatures), row.names(permhist.wide))){
      stop("SVC permutation features do not match the common gene universe")
    }
    permhist.wide = cbind(
      permhist.wide,
      "SVC RFE-CV candidates" =
        SVC.permfeatures[row.names(permhist.wide)]
    )
  }

  if("$params.enetstep" == "true"){
    enet.permfeatures = read_permfeatures_enet(
      "enet",
      row.names(permhist.wide)
    )
    permhist.wide = cbind(
      permhist.wide,
      "Stability-selected elastic net" =
        enet.permfeatures[row.names(permhist.wide)]
    )
  }

  permhist.in = melt(permhist.wide)

  # Generate breaks contingent on number of perms
  cutbreaks = cut(seq(1, nperms + 1), breaks = 5, right = TRUE) %>% levels() %>%
    str_match("\\\\([^,]*") %>%
    str_remove("\\\\(") %>%
    as.numeric() %>%
    round()
  cutbreaks[1] = 1

  gg.permhist = ggplot(permhist.in, aes(x = value)) +
    geom_histogram(color = "black", fill = "white", bins = 100) +
    scale_x_continuous(limits = c(-1, (nperms + 1)), breaks = cutbreaks) +
    labs(x = "Number of permutations in which each gene was selected",
         y = "Number of genes") +
    theme_bw() +
    facet_grid(~Var2) +
    theme(strip.background = element_blank())

  ggsave(gg.permhist,
   filename = "permute_hist.pdf",
   device = "pdf", bg = "transparent",
   width = 30, height = 20, units = "cm")

  # Compare model fit for rarely and frequently selected genes
  if(nperms >= 5){

    floor   = nperms / 1000
    ceiling = nperms / 5

    permhist.deseq  = subset(permhist.in, Var2 == "DESeq2")
    deseq.trust     = subset(permhist.deseq, value <= floor)
    deseq.untrust   = subset(permhist.deseq, value >= ceiling)

    deseq.poorfits = get(load("$deseq_poorfits"))

    # NB dummy NAs prevent errors when there are 0 (un)trustworthy genes
    permhist.edger = subset(permhist.in, Var2 == "edgeR")
    edger.trust    = subset(permhist.edger, value <= floor)
    edger.untrust  = subset(permhist.edger, value >= ceiling)

    trustcomp = rbind(
      data.frame(poorfit = c(deseq.poorfits[deseq.trust\$Var1],   NA), val = "Trustworthy",   method = "DESeq2"),
      data.frame(poorfit = c(deseq.poorfits[deseq.untrust\$Var1], NA), val = "Untrustworthy", method = "DESeq2"),
      data.frame(poorfit = c(deseq.poorfits[edger.trust\$Var1],   NA), val = "Trustworthy",   method = "edgeR"),
      data.frame(poorfit = c(deseq.poorfits[edger.untrust\$Var1], NA), val = "Untrustworthy", method = "edgeR")
    )

    gg.deseq.trust =
      ggplot(trustcomp, aes(x = as.factor(val), y = poorfit)) +
      geom_violin() +
      labs(x = "Gene reliability",
           y = "Poorness of fit to Negative Binomial distribution") +
      stat_compare_means(method = "wilcox.test",
                         label = "p.signif",
                         size = 5,
                         comparisons = list(c("Trustworthy", "Untrustworthy"))) +
      facet_wrap(~ method) +
      theme_bw()

    ggsave(gg.deseq.trust,
          filename = "permute_violins.pdf",
          device = "pdf", bg = "transparent",
          width = 20, height = 12, units = "cm")

  } else {
    message(
      "Skipping permute_violins.pdf: this diagnostic requires at least five permutations; ",
      "this run used ", as.integer(nperms), "."
    )
  }
  """
}
