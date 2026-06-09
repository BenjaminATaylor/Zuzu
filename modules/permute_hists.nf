process PERMUTE_HISTS{

  publishDir "$params.outdir"

  input:
  path(deseq_infiles,  stageAs: 'deseq/deseq_table_?.csv')
  path(edger_infiles,  stageAs: 'edger/edger_table_?.csv')
  path(wilcox_infiles, stageAs: 'wilcox/wilcox_table_?.csv')
  path(SVC_infiles,    stageAs: 'svc/svc_table_?.csv')
  path deseq_poorfits

  output:
  path 'permute_hist.pdf'
  path 'permute_violins.pdf'

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse", quietly = TRUE)
  library("reshape2", quietly = TRUE)
  library("ggpubr")

  nperms = as.numeric("$params.nperms")

  read_permdegs_r = function(dir){
    lapply(list.files(dir, full.names = TRUE), function(x){
      res <- read.csv(x, row.names = 1)
      data.frame(result = (res\$padj < 0.05), row.names = row.names(res))
    }) %>% list_cbind() %>% rowSums(na.rm = TRUE)
  }

  read_permdegs_py = function(dir){
    lapply(list.files(dir, full.names = TRUE), function(x){
      res <- read.csv(x, row.names = 1)
      data.frame(result = (res\$DEGstatus == 1), row.names = row.names(res))
    }) %>% list_cbind() %>% rowSums(na.rm = TRUE)
  }

  deseq.permdegs  = read_permdegs_r("deseq")
  edger.permdegs  = read_permdegs_r("edger")
  wilcox.permdegs = read_permdegs_r("wilcox")

  if("$params.mlstep" == "true"){
    SVC.permdegs = read_permdegs_py("svc")
    permhist.in = cbind(deseq.permdegs,
                        edger.permdegs,
                        wilcox.permdegs,
                        SVC.permdegs) %>%
      `colnames<-`(c("DESeq2", "edgeR", "Wilcoxon", "SVC")) %>%
      melt()
  } else {
    permhist.in = cbind(deseq.permdegs,
                        edger.permdegs,
                        wilcox.permdegs) %>%
      `colnames<-`(c("DESeq2", "edgeR", "Wilcoxon")) %>%
      melt()
  }

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
    labs(x = "Number of permutations in which gene was\nincorrectly identified as a DEG",
         y = "Number of genes") +
    theme_bw() +
    facet_grid(~Var2) +
    theme(strip.background = element_blank())

  ggsave(gg.permhist,
   filename = "permute_hist.pdf",
   device = "pdf", bg = "transparent",
   width = 30, height = 20, units = "cm")

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
    dummy = NA
    save(dummy, file = "permute_violins.pdf")
  }
  """
}
