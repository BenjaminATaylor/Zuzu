process PERMUTE_HISTS{

  publishDir "$params.outdir"

  input: 
  val deseq_infiles
  val edger_infiles
  val wilcox_infiles
  val SVC_infiles
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

  # Function to read in tables produced by R
  readfun = function(x){
    res <- read.csv(x,row.names = 1)
    df <- data.frame(result = (res\$padj<0.05), row.names = row.names(res))
    return(df)
  }

  # Function to read in files produced by Python
  pythout_readfun = function(x){
    res <- read.csv(x,row.names = 1)
    df <- data.frame(result = (res\$DEGstatus==1), row.names = row.names(res))
    return(df)
  }

  deseq.inlist = str_remove_all("$deseq_infiles","[\\\\[\\\\] ]") %>% 
    strsplit(split = ",") %>% unlist()
  deseq.permdegs = lapply(deseq.inlist,readfun) %>% 
    list_cbind() %>%
    rowSums(na.rm=TRUE)

  edger.inlist = str_remove_all("$edger_infiles","[\\\\[\\\\] ]") %>% 
    strsplit(split = ",") %>% unlist()
  edger.permdegs = lapply(edger.inlist,readfun) %>% 
    list_cbind() %>%
    rowSums(na.rm=TRUE)

  wilcox.inlist = str_remove_all("$wilcox_infiles","[\\\\[\\\\] ]") %>% 
    strsplit(split = ",") %>% unlist()
  wilcox.permdegs = lapply(wilcox.inlist,readfun) %>% 
    list_cbind() %>%
    rowSums(na.rm=TRUE)

  # Only include ML outputs if specified by user
  if("$params.mlstep" == "true"){
    SVC.inlist = str_remove_all("$SVC_infiles","[\\\\[\\\\] ]") %>% 
      strsplit(split = ",") %>% unlist()
    SVC.permdegs = lapply(SVC.inlist,pythout_readfun) %>% 
      list_cbind() %>%
      rowSums(na.rm=TRUE)
    
    permhist.in = cbind(deseq.permdegs,
                        edger.permdegs,
                        wilcox.permdegs,
                        SVC.permdegs) %>%
      `colnames<-`(c("DESeq2",
                    "edgeR",
                    "Wilcoxon",
                    "SVC")) %>%
      melt()
  } else {
    permhist.in = cbind(deseq.permdegs,
                        edger.permdegs,
                        wilcox.permdegs) %>%
      `colnames<-`(c("DESeq2",
                    "edgeR",
                    "Wilcoxon")) %>%
      melt()
}


  # Generate breaks contingent on number of reps
  cutbreaks = cut(seq(1,nperms+1),breaks = 5, right = TRUE) %>% levels() %>%
    str_match("\\\\([^,]*") %>%
    str_remove("\\\\(") %>%
    as.numeric() %>%
    round()
  #deal with an idiosyncrasy that causes cut to return 0 as the start value of the first bin instead of 1
  cutbreaks[1] = 1
  
  #plot
  gg.permhist = ggplot(permhist.in, aes(x = value)) +
    geom_histogram(color = "black", fill = "white", bins = 100) +
    scale_x_continuous(limits = c(-1,(nperms+1)), breaks = cutbreaks) +
    labs(x = "Number of permutations in which gene was\nincorrectly identified as a DEG",
         y = "Number of genes") +
    theme_bw() +
    facet_grid(~Var2) +
    theme(strip.background = element_blank())

  ggsave(gg.permhist, 
   filename = "permute_hist.pdf",
   device = "pdf", bg = "transparent",
   width =  30, height = 20, units = "cm")
   
  # Also here we can run the check where we see if frequently misidentified genes tend to be nonbinomial
  # This is only really meaningful when we have a large-ish number of permutations, but we'll do it with fewer for testing purposes
  if(nperms >= 5){
    
    floor = nperms/1000 # Genes that are misattributed less than this often are 'trustworthy'
    ceiling = nperms/5 # Genes that are misattributed more than this often are 'untrustworthy'
    
    permhist.deseq = subset(permhist.in, Var2 == "DESeq2")
    deseq.trust = subset(permhist.deseq, value<=floor)
    deseq.untrust = subset(permhist.deseq, value>=ceiling)
    
    #load poorness-of-fit data for deseq2 
    deseq.poorfits = get(load("$deseq_poorfits"))
    
    # NB dummy NAs are appended here to prevent errors occurring when there are 0 '(un)trustworthy' genes
    deseq.trust = rbind(
      data.frame(poorfit = c(deseq.poorfits[deseq.trust\$Var1],NA), val = "Trustworthy"),
      data.frame(poorfit = c(deseq.poorfits[deseq.untrust\$Var1],NA), val = "Untrustworthy")
    ) %>% mutate(method = "DESeq2")
    
    permhist.edger = subset(permhist.in, Var2 == "edgeR")
    edger.trust = subset(permhist.edger, value<=floor)
    edger.untrust = subset(permhist.edger, value>=ceiling)
    
    #TODO: load poorness-of-fit data for edger instead of using DESeq fits here
    #edger.poorfits = get(load("/Users/benjamin/Repositories/Zuzu/work/c5/cbd7b1c7e36c0109c2bff8a91b63be/edger_poorfits.RData"))
    
    edger.trust = rbind(
      data.frame(poorfit = c(deseq.poorfits[edger.trust\$Var1],NA), val = "Trustworthy"),
      data.frame(poorfit = c(deseq.poorfits[edger.untrust\$Var1],NA), val = "Untrustworthy")
    ) %>% mutate(method = "edgeR")
    
    trustcomp = rbind(deseq.trust,edger.trust)
    
    gg.deseq.trust = 
      ggplot(trustcomp, aes(x = as.factor(val), y = poorfit)) +
      geom_violin() + 
      labs(x = "Gene reliability",
          y = "Poorness of fit to Negative Binomial distribution") +
      stat_compare_means(method = "wilcox.test",
                        label = "p.signif",
                        size = 5,
                        comparisons = list(c("Trustworthy","Untrustworthy"))) + #,
    #label.y = c(0.31, 0.275, 0.255) +
    facet_wrap(~ method) +
    theme_bw() 

    ggsave(gg.deseq.trust, 
          filename = "permute_violins.pdf",
          device = "pdf", bg = "transparent",
          width =  20, height = 12, units = "cm")

  } else {
  # If not enough perms, fool nextflow with dummy output 
    dummy = NA
    save(dummy , file = "permute_violins.pdf")
  }

  """
}