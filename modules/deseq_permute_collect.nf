process DESEQ_PERMUTE_COLLECT{

  //debug true

  input: 
  val nDEGs_list

  output:
  path "deseq_nDEGs_list.RData"

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse")

  objlist = str_remove_all("$nDEGs_list","[\\\\[\\\\] ]") %>% strsplit(split = ",") %>% unlist()
  permdeglist = sapply(objlist, function(x) get(load(x))) %>% unname()

  save(permdeglist, file="deseq_nDEGs_list.RData")
  """
}