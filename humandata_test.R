library(tidyverse)
library(DESeq2)

# Import real testing dataset
humancounts = read.csv("input/test3/GSE91061_BMS038109Sample.hg19KnownGene.raw.csv", row.names = 1, check.names = T)
humandata = read.csv("input/test3/data_phenotype.csv") %>% select(c("sample","phenotype")) %>%
  mutate(sample = str_replace(sample, "-", "."))

# subset to just 60 samples and 1000 genes
groupsize = 30
picksamples = c(subset(humandata, phenotype == "Pre-treatment")$sample[1:groupsize],
                subset(humandata, phenotype == "On-treatment")$sample[1:groupsize])

#test keeping all samples
picksamples = c(humandata$sample)

genesize = nrow(humancounts)
pickgenes = rownames(humancounts)[1:genesize]
testcounts = humancounts[pickgenes,picksamples]
testdata = subset(humandata, sample %in% picksamples)
#check conformity of new dataframes
table(testdata$sample == colnames(testcounts))

# basic DESeq2 design
dds.gene = DESeqDataSetFromMatrix(countData = humancounts,
                                  colData = humandata,
                                  design = as.formula(~phenotype))
dds.gene.deg = DESeq(dds.gene)

#check data look reasonable
boxplot(log10(assays(dds.gene.deg)[["cooks"]]), range=0, las=2)
boxplot(log10(assays(dds.gene.deg)[["counts"]]), range=0, las=2)
plotDispEsts(dds.gene.deg)

# Okay that seems in the right ballpark for number of DEGs
table(results(dds.gene.deg)$padj<0.05)
thesedegs = row.names(subset(results(dds.gene.deg),padj<0.05))


BiocManager::install("seqgendiff", Ncpus = 6)


library(airway)
data("airway")
library(seqgendiff)

coldat <- colData(airway)[, c("cell", "dex")]
allzero <- rowSums(assay(airway)) < 10^-6
airway <- airway[!allzero, ]

thout <- thin_2group(mat = assay(airway), 
                     prop_null = 0.9, 
                     signal_fun = stats::rnorm,
                     signal_params = list(mean = 0, sd = 0.8))

X <- cbind(thout$design_obs, thout$designmat)
dim(thout$mat)

