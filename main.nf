#!/usr/bin/env nextflow

// Modules to include
include { CLEANINPUTS } from './modules/cleaninputs.nf'
include { QUALITYCONTROL } from './modules/qualitycontrol.nf'

//exit 1, 'DEBUG'

// Validate inputs
if (params.samplesheet) { ch_samplesheet = file(params.samplesheet) } else { exit 1, 'Input samplesheet not specified!' }
if (params.countsframe) { ch_countsframe = file(params.countsframe) } else { exit 1, 'Input count frame not specified!' }
if (params.reflevel) { ch_reflevel = params.reflevel } else { exit 1, 'Reference level not specified!' }

// Set number of permutations for fakey datasets
params.nperms = 7

println("Sample sheet: " + ch_samplesheet)
println("Counts matrix: " + ch_countsframe)
println("Reference level: " + ch_reflevel)

process DESEQ_BASIC{

  input: 
  tuple path(samplesheet), path(countsframe) 

  output:
  path "dds_deg.RData"

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse")
  library("DESeq2")

  samplesheet = read.csv("$samplesheet")
  countsframe.clean = read.csv("$countsframe", row.names = 1)

  #generate model
  dds = DESeqDataSetFromMatrix(countData = countsframe.clean,
                               colData = samplesheet,
                               design = as.formula(~phenotype))
  # Run the default analysis for DESeq2
  dds.deg = DESeq(dds, fitType = "parametric", betaPrior = FALSE)
  save(dds.deg, file = "dds_deg.RData")
  """
}

process DATA_PERMUTE{

  input: 
  tuple path(samplesheet), path(countsframe)
  val x

  output:
  tuple path("samplesheet_perm.csv"), path("countsframe_perm.csv")

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse")

  #print(paste0("DEBUG: ", "$samplesheet", " ", "$x"))

  samplesheet = read.csv("$samplesheet")
  countsframe.clean = read.csv("$countsframe", row.names = 1)

  # Permute colnames at random while preserving gene counts
  set.seed($x)
  permcols = sample(colnames(countsframe.clean),replace = FALSE)
  countsframe.perm = countsframe.clean %>% `colnames<-`(permcols)
  samplesheet.perm = samplesheet[match(samplesheet\$sample, permcols),]

  write.csv(samplesheet.perm,"samplesheet_perm.csv", row.names = FALSE)
  write.csv(countsframe.perm,"countsframe_perm.csv")

  """
}

process DESEQ_PERMUTE{

  input: 
  tuple path(samplesheet), path(countsframe) 

  output:
  path "dds_deg.RData"
  path "nDEGs.RData", emit: nDEGs


  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse")
  library("DESeq2")

  samplesheet = read.csv("$samplesheet")
  countsframe.clean = read.csv("$countsframe", row.names = 1)

  #generate model
  dds = DESeqDataSetFromMatrix(countData = countsframe.clean,
                               colData = samplesheet,
                               design = as.formula(~phenotype))
  # Run the default analysis for DESeq2
  dds.deg = DESeq(dds, fitType = "parametric", betaPrior = FALSE)
  nDEGs = nrow(subset(results(dds.deg),padj<0.05))
  save(nDEGs,file = "nDEGs.RData")
  save(dds.deg, file = "dds_deg.RData")
  """
}

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





workflow {
  CLEANINPUTS(ch_samplesheet, ch_countsframe, ch_reflevel)
  QUALITYCONTROL(CLEANINPUTS.out)
  DESEQ_BASIC(CLEANINPUTS.out)

  perms = Channel.from(1..params.nperms)
  DATA_PERMUTE(CLEANINPUTS.out, perms) |
  DESEQ_PERMUTE 

  DESEQ_PERMUTE.out.nDEGs.collect() |
  DESEQ_PERMUTE_COLLECT 
  
  DESEQ_PERMUTE_COLLECT.out.view()
}
