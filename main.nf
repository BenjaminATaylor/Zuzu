#!/usr/bin/env nextflow

// Modules to include
include { CLEANINPUTS } from './modules/cleaninputs.nf'
//exit 1, 'DEBUG'

// Validate inputs
if (params.samplesheet) { ch_samplesheet = file(params.samplesheet) } else { exit 1, 'Input samplesheet not specified!' }
if (params.countsframe) { ch_countsframe = file(params.countsframe) } else { exit 1, 'Input count frame not specified!' }
if (params.reflevel) { ch_reflevel = params.reflevel } else { exit 1, 'Reference level not specified!' }

println("Sample sheet: " + ch_samplesheet)
println("Counts matrix: " + ch_countsframe)
println("Reference level: " + ch_reflevel)

process QUALITYCONTROL {

  input: 
  tuple path(samplesheet), path(countsframe)

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse")



  """

}

workflow {
  CLEANINPUTS(ch_samplesheet, ch_countsframe, ch_reflevel)
    .view()
}
