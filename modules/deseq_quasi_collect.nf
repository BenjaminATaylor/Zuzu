process DESEQ_QUASI_COLLECT{

  debug true

  input: 
  val 'quasiout'

  output:
  tuple path("quasi_outframe.csv"), val("DESeq2")

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse")

  # print("$quasiout")

  objlist = str_remove_all("$quasiout","[\\\\[\\\\] ]") %>% strsplit(split = ",") %>% unlist()
  quasi.outframe = sapply(objlist, function(x) read.csv(x)) %>% 
    t() %>% `rownames<-`(NULL) %>% 
    data.frame() %>% mutate_all(as.numeric)

  write.csv(quasi.outframe, file = "quasi_outframe.csv", row.names = FALSE)
  """
}