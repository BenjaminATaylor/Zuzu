process DESEQ_PERMUTE_COLLECT{

  //debug true

  input: 
  val nDEGs_list

  output:
  tuple path("nDEGs_list.RData"), val("DESeq2")

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse")

  objlist = str_remove_all("$nDEGs_list","[\\\\[\\\\] ]") %>% strsplit(split = ",") %>% unlist()
  permdeglist = sapply(objlist, function(x) get(load(x))) %>% unname()

  save(permdeglist, file="nDEGs_list.RData")
  """
}