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

process DESEQ_ORIGINAL{

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

process DESEQ_PERMUTE{

  input: 
  tuple path(samplesheet), path(countsframe)
  val x

  output:
  path "dds_deg_permute.RData"

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse")
  library("DESeq2")

  #print(paste0("DEBUG: ", "$samplesheet", " ", "$x"))

  samplesheet = read.csv("$samplesheet")
  countsframe.clean = read.csv("$countsframe", row.names = 1)

  # Permute colnames at random while preserving gene counts
  set.seed($x)
  permcols = sample(colnames(countsframe.clean),replace = FALSE)
  countsframe.perm = countsframe.clean %>% `colnames<-`(permcols)
  samplesheet.perm = samplesheet[match(samplesheet\$sample, permcols),]

  #generate model
  dds.perm = DESeqDataSetFromMatrix(countData = countsframe.perm,
                                    colData = samplesheet.perm,
                                    design = as.formula(~phenotype))
  # Run the default analysis for DESeq2
  dds.deg.perm = DESeq(dds.perm, fitType = "parametric", betaPrior = FALSE)

  # Save the *number* of DEGs as output (that's all we care about for now, I think?)
  save(dds.deg.perm, file = "dds_deg_permute.RData")
  """
}



workflow {
  CLEANINPUTS(ch_samplesheet, ch_countsframe, ch_reflevel)
  QUALITYCONTROL(CLEANINPUTS.out)
  DESEQ_ORIGINAL(CLEANINPUTS.out)

  perms = Channel.from(1..params.nperms)
  DESEQ_PERMUTE(CLEANINPUTS.out, perms)

  DESEQ_PERMUTE.out.view()
}
