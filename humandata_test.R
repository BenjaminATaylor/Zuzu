library(tidyverse)
library(DESeq2)

# Import real testing dataset
humancounts = read.csv("input/test3/GSE91061_BMS038109Sample.hg19KnownGene.raw.csv", row.names = 1, check.names = F)
humandata = read.csv("input/test3/data_phenotype.csv") %>% select(c("sample","phenotype"))

# subset to just 60 samples and 1000 genes
groupsize = 30
picksamples = c(subset(humandata, phenotype == "Pre-treatment")$sample[1:groupsize],
                subset(humandata, phenotype == "On-treatment")$sample[1:groupsize])
genesize = 1000
pickgenes = rownames(humancounts)[1:genesize]
testcounts = humancounts[pickgenes,picksamples]
testdata = subset(humandata, sample %in% picksamples)
#check conformity of new dataframes
table(testdata$sample == colnames(testcounts))

