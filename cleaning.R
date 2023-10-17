##### Data import and cleaning #####

# Libs
library(tidyverse)

# Setup parameters
setup.workdir = "/Users/benjamin/Repositories/Zuzu/"
setup.samplesheet = "input/test/samplesheet_phenotype.csv"
setup.countsframe = "input/test/salmon.merged.gene_counts.tsv"
setup.reflevel = "nurse"

# Setup
setwd(setup.workdir)
samplesheet = read.csv(setup.samplesheet)
countsframe = read.delim(setup.countsframe, row.names = 1) %>% select_if(is.numeric)

# Check count and phenotype matrix conformity
stopifnot(identical(sort(samplesheet$sample), sort(colnames(countsframe))))
# Order both frames
samplesheet = samplesheet[order(samplesheet$sample),]
countsframe = countsframe[,order(colnames(countsframe))]

# Check that the 'phenotype' conditions has exactly two conditions, and set a reference level
reflevel = setup.reflevel
stopifnot(reflevel %in% levels(as.factor(samplesheet$phenotype)))
samplesheet$phenotype = as.factor(samplesheet$phenotype) %>% relevel(ref = reflevel)
stopifnot(length(levels(samplesheet$phenotype))==2)
