process PERMUTE_PLOTS{

  debug false

  publishDir "$params.outdir"

  input:
  path deseq_table
  val deseq_perms
  path edger_table
  val edger_perms
  path wilcox_table
  val wilcox_perms
  path SVC_table
  val SVC_perms

  output:
  path "permute_plot.pdf"

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse", quietly = TRUE)
  library("DESeq2", quietly = TRUE)

  ## DESeq inputs
  # DEGs from full model
  deseq.table = read.csv("$deseq_table", row.names = 1)
  deseq.ndegs = nrow(subset(deseq.table, padj<0.05))
  # Permutations
  deseq.inlist = str_remove_all("$deseq_perms","[\\\\[\\\\] ]") %>% 
    strsplit(split = ",") %>% unlist()
  deseq.perms = sapply(deseq.inlist, function(x) get(load(x))) %>% unname()
  # Collated input
  deseq.permute.input = data.frame( method = "DESeq2", 
                                  nSelected = deseq.ndegs,
                                  permutes = deseq.perms
                                  )

  ## edgeR inputs
  # DEGs from full model
  edger.table = read.csv("$edger_table", row.names = 1)
  edger.ndegs = nrow(subset(edger.table, padj<0.05))
  # Permutations
  edger.inlist = str_remove_all("$edger_perms","[\\\\[\\\\] ]") %>% 
    strsplit(split = ",") %>% unlist()
  edger.perms = sapply(edger.inlist, function(x) get(load(x))) %>% unname()
  # Collated input
  edger.permute.input =  data.frame(method = "edgeR", 
                                    nSelected = edger.ndegs,
                                    permutes = edger.perms)

  ## wilcoxon inputs
  # DEGs from full model
  wilcox.table = read.csv("$wilcox_table", row.names = 1)
  wilcox.ndegs = nrow(subset(wilcox.table, padj<0.05))
  # Permutations
  wilcox.inlist = str_remove_all("$wilcox_perms","[\\\\[\\\\] ]") %>% 
    strsplit(split = ",") %>% unlist()
  wilcox.perms = sapply(wilcox.inlist, function(x) get(load(x))) %>% unname()
  # Collated input
  wilcox.permute.input =  data.frame(method = "wilcox", 
                                    nSelected = wilcox.ndegs,
                                    permutes = wilcox.perms)
  
  # Only include ML outputs if specified by user
  if("$params.mlstep" == "true"){
    ## SVC inputs
    # Candidate features from the full model
    SVC.table = read.csv("$SVC_table", row.names = 1)
    SVC.nfeatures = nrow(subset(SVC.table, selected_feature == 1))
    # Permutations
    SVC.inlist = str_remove_all("$SVC_perms","[\\\\[\\\\] ]") %>% 
      strsplit(split = ",") %>% unlist()
    SVC.perms = sapply(SVC.inlist, function(x) as.integer(read.delim(x, header = FALSE))) %>% unname()
    # Collated input
    SVC.permute.input =  data.frame(method = "SVC RFE-CV candidates",
                                    nSelected = SVC.nfeatures,
                                    permutes = SVC.perms)

    # Combined inputs for plotting
    permuteplot.input = rbind(deseq.permute.input,
                              edger.permute.input,
                              wilcox.permute.input,
                              SVC.permute.input)
  } else {
    # Combined inputs for plotting
    permuteplot.input = rbind(deseq.permute.input,
                              edger.permute.input,
                              wilcox.permute.input)
    
  }

  print(gg.permute <- ggplot(permuteplot.input, aes(x = method, y = permutes)) +
          geom_point(size = 3, alpha = 0.7) +
          geom_point(aes(x = method, y = nSelected),
                    size =4 , color = "black", fill = "red", shape = 23) +
          labs(x = "Method", y = "Number of selected genes") +
          theme_bw()
  )

  ggsave(gg.permute, 
       filename = "permute_plot.pdf",
       device = "pdf", bg = "transparent",
       width =  30, height = 20, units = "cm")

  """

}
