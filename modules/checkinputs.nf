process CHECKINPUTS {

    conda 'r-tidyverse-1.3.1'

    input: 
    path samplesheet
    path countsframe
    val reflevel

    output:
    tuple path("samplesheet.csv"), path("countsframe.csv")

    script:
    """
    #!/usr/bin/env Rscript
    library("tidyverse")

    # Setup
    samplesheet = read.csv("$samplesheet")
    countsframe = read.delim("$countsframe", row.names = 1) %>% select_if(is.numeric)

    # Check count and phenotype matrix conformity
    stopifnot(identical(sort(samplesheet\$sample), sort(colnames(countsframe))))
    # Order both frames
    samplesheet = samplesheet[order(samplesheet\$sample),]
    countsframe = countsframe[,order(colnames(countsframe))]

    # Check that the 'phenotype' conditions has exactly two conditions, and set a reference level
    stopifnot("$reflevel" %in% levels(as.factor(samplesheet\$phenotype)))
    samplesheet\$phenotype = as.factor(samplesheet\$phenotype) %>% relevel(ref = "$reflevel")
    stopifnot(length(levels(samplesheet\$phenotype))==2)
    
    write.csv(samplesheet,"samplesheet.csv", row.names = FALSE)
    write.csv(countsframe,"countsframe.csv")
    """
}