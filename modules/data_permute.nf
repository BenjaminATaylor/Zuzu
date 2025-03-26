process DATA_PERMUTE{

  label 'small_job'

  input: 
  tuple path(samplesheet), path(countsframe)
  val x

  output:
  tuple path(samplesheet), path("countsframe_perm.csv"), emit: frames
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
  permute_fun = function(x){
    for(i in 1:100){ # NB if we don't get a proper permutation after 100 runs, we give up because the sample sizes are obviously inappropriate
      # Permute and check whether any of the newly permuted rows have perfectly recaptured the original groupings
      permrow = sample(x,replace = FALSE) 
      newcols = samplesheet\$phenotype[match(names(permrow), samplesheet\$sample)]
      proptrue = sum(newcols == samplesheet\$phenotype)/length(newcols)
      # If more than 80% of sample labels have remained the same or have flipped, re-permute everything
      if(abs(proptrue-0.5)<0.3){
        break
      }else if(i==100){
        stop("Unable to achieve a balanced permutation. Are your sample sizes large enough?")
      }
    }
    return(data.frame(t(permrow)))
  }

  #coerce back to df
  countsframe.perm = apply(countsframe.clean, 1, permute_fun) %>% 
    bind_rows(., .id = "column_label") %>%
    column_to_rownames(var = "column_label") %>%
    `colnames<-`(colnames(countsframe.clean))

  #write.csv(samplesheet,"samplesheet_perm.csv", row.names = FALSE)
  write.csv(countsframe.perm,"countsframe_perm.csv")

  """
}