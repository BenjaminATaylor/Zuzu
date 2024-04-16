process CLEANINPUTS {

    label 'small_job'

    input: 
    path samplesheet
    path countsframe
    val reflevel

    output:
    tuple path("samplesheet.csv"), path("countsframe_clean.csv")

    script:
    """
    #!/usr/bin/env Rscript
    library("tidyverse", quietly = TRUE)

    # Setup
    samplesheet = read.csv("$samplesheet")
    countsframe = read.delim("$countsframe", row.names = 1, check.names = FALSE) %>% 
        select_if(is.numeric)

    # Check count and phenotype matrix conformity
    stopifnot(identical(sort(samplesheet\$sample), sort(colnames(countsframe))))
    # Order both frames
    samplesheet = samplesheet[order(samplesheet\$sample),]
    countsframe = countsframe[,order(colnames(countsframe))]

    # Check that the 'phenotype' conditions has exactly two conditions, and set a reference level
    stopifnot("$reflevel" %in% levels(as.factor(samplesheet\$phenotype)))
    samplesheet\$phenotype = as.factor(samplesheet\$phenotype) %>% relevel(ref = "$reflevel")
    stopifnot(length(levels(samplesheet\$phenotype))==2)

    # Round all values to integers
    countsframe = round(countsframe)

    # Now clean for low counts. Require mean per-sample counts greater than 5 in at least one group
    checkmeans = function(x){
    thesesamples = subset(samplesheet, phenotype == x)\$sample
    rowMeans(countsframe[,thesesamples])>5
    }
    keeprows = row.names(countsframe)[rowSums(sapply(levels(samplesheet\$phenotype),checkmeans))>0]
    countsframe.clean = countsframe[keeprows,]

    write.csv(samplesheet,"samplesheet.csv", row.names = FALSE)
    write.csv(countsframe.clean,"countsframe_clean.csv")
    """
}