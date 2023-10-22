#!/usr/bin/env nextflow

// Modules to include
include { CLEANINPUTS } from './modules/cleaninputs.nf'
include { QUALITYCONTROL } from './modules/qualitycontrol.nf'
include { DESEQ_BASIC } from './modules/deseq_basic.nf'
include { DATA_PERMUTE } from './modules/data_permute.nf'
include { DESEQ_PERMUTE } from './modules/deseq_permute.nf'
include { DESEQ_PERMUTE_COLLECT } from './modules/deseq_permute_collect.nf'
include { DATA_QUASI } from './modules/data_quasi.nf'
include { DESEQ_QUASI } from './modules/deseq_quasi.nf'
include { DESEQ_QUASI_COLLECT } from './modules/deseq_quasi_collect.nf'

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

process EDGER_BASIC{

  input: 
  tuple path(samplesheet), path(countsframe) 

  output:
  path "edger_table.csv"

  script:
  """
  #!/usr/bin/env Rscript
  library("edgeR") 

  # Input data
  samplesheet = read.csv("$samplesheet")
  countsframe.clean = read.csv("$countsframe", row.names = 1)

  # Basic edgeR workflow, taken from the documentation
  dds.edge = DGEList(counts=countsframe.clean,group=samplesheet\$phenotype)
  dds.edge = normLibSizes(dds.edge)
  design.edge = model.matrix(~samplesheet\$phenotype)
  dds.edge.est = estimateDisp(dds.edge,design.edge)
  edge.fit = glmQLFit(dds.edge.est,design.edge)
  edge.qlf = glmQLFTest(edge.fit,coef=2)

  # Tabular output (rather than topTags, which only gives the most highly-DE results)
  edge.out = edge.qlf\$table
  write.csv(edge.out, row.names = T, file="edger_table.csv")
  """
}


workflow {
  //Initial data QC and cleanup
  CLEANINPUTS(ch_samplesheet, ch_countsframe, ch_reflevel)
  QUALITYCONTROL(CLEANINPUTS.out)

  //Basic analysis for each method
  DESEQ_BASIC(CLEANINPUTS.out)
  EDGER_BASIC(CLEANINPUTS.out)

  //Permutation analysis with no true signal
  perms = Channel.from(1..params.nperms)
  DATA_PERMUTE(CLEANINPUTS.out, perms)
  DESEQ_PERMUTE(DATA_PERMUTE.out)
  DESEQ_PERMUTE.out.nDEGs.view()


  DESEQ_PERMUTE.out.nDEGs
  .collect() |
  DESEQ_PERMUTE_COLLECT 
  
  DATA_QUASI(DESEQ_BASIC.out, CLEANINPUTS.out, perms) |
  DESEQ_QUASI 

  DESEQ_QUASI.out
  .collect() |
  DESEQ_QUASI_COLLECT

  DESEQ_QUASI_COLLECT.out.view()

}
