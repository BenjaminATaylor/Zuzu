process CREATE_BREAKS {

  debug true

  input:
  path samplesheet

  output:
  file 'breaks.csv'
  //env refouts, emit: refouts
  //env altouts, emit: altouts

  """
  #!/usr/bin/env Rscript
  library("tidyverse", quietly = TRUE)
  samplesheet = read.csv("$samplesheet")
  reflevel = "$params.reflevel"

  reftot = nrow(subset(samplesheet, phenotype == reflevel))
  alttot = nrow(subset(samplesheet, phenotype != reflevel))

  refbreaks = cut(seq(1,reftot),breaks = 3, right = TRUE) %>% levels() %>%
    str_match("\\\\,[^]]*") %>%
    str_remove(",") %>%
    as.numeric() %>%
    round()

  altbreaks = cut(seq(1,alttot),breaks = 3, right = TRUE) %>% levels() %>%
    str_match("\\\\,[^]]*") %>%
    str_remove(",") %>%
    as.numeric() %>%
    round()

  outbreaks = data.frame(refbreaks, altbreaks)
  write.table(outbreaks, "breaks.csv", sep = ",", row.names = F, col.names = F)
  """
}