#!/usr/bin/env nextflow

// Modules to include
//include { CREATESAMPLESHEET } from './modules/createsamplesheet.nf'
//exit 1, 'DEBUG'

// Validate inputs
if (params.samplesheet) { ch_samplesheet = file(params.samplesheet) } else { exit 1, 'Input samplesheet not specified!' }
if (params.countsframe) { ch_countsframe = file(params.countsframe) } else { exit 1, 'Input count frame not specified!' }
if (params.reflevel) { ch_reflevel = params.reflevel } else { exit 1, 'Reference level not specified!' }

println("Sample sheet: " + ch_samplesheet)
println("Counts matrix: " + ch_countsframe)
println("Reference level: " + ch_reflevel)

process CHECKDATA {

    conda 'r-tidyverse-1.3.1'

    input: 
    path samplesheet
    path countsframe
    val reflevel

    output:
    stdout

    script:
    """
    #!/usr/bin/env Rscript
    library("tidyverse")

    # Setup
    samplesheet = read.csv("$samplesheet")
    countsframe = read.delim("$countsframe", row.names = 1) %>% select_if(is.numeric)

    nrow(countsframe)
    """
}

workflow {
  CHECKDATA(ch_samplesheet, ch_countsframe, ch_reflevel)
    .view()
}
