process DATA_PERMUTE{

  input: 
  tuple path(samplesheet), path(countsframe)
  val x

  output:
  tuple path("samplesheet_perm.csv"), path("countsframe_perm.csv"), emit: frames
  val x

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse")

  #print(paste0("DEBUG: ", "$samplesheet", " ", "$x"))

  samplesheet = read.csv("$samplesheet")
  countsframe.clean = read.csv("$countsframe", row.names = 1, check.names = FALSE)

  # Permute colnames at random while preserving gene counts
  set.seed($x)
  permcols = sample(colnames(countsframe.clean),replace = FALSE)
  countsframe.perm = countsframe.clean %>% `colnames<-`(permcols)
  samplesheet.perm = samplesheet[match(samplesheet\$sample, permcols),]

  write.csv(samplesheet.perm,"samplesheet_perm.csv", row.names = FALSE)
  write.csv(countsframe.perm,"countsframe_perm.csv")

  """
}