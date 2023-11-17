library(compcodeR)
library(tidyverse)
library(DESeq2)

setwd("/Users/benjamin/Repositories/Zuzu/")

# generate synthetic data with 10000 genes, 20 samples/group, and 1200 differentially expressed genes
synthtest = generateSyntheticData(dataset = "testset_1", n.vars = 5000, 
                                  samples.per.cond = 50, n.diffexp = 2500, 
                                  repl.id = 1, seqdepth = 1e6, 
                                  fraction.upregulated = 0.5, 
                                  between.group.diffdisp = FALSE, 
                                  filter.threshold.total = 10, 
                                  filter.threshold.mediancpm = 0, 
                                  fraction.non.overdispersed = 0, 
                                  output.file = NULL)

# clean data
synthtest.counts = synthtest@count.matrix
synthtest.metadata = 
  synthtest@sample.annotations %>% 
  mutate(phenotype = as.factor(condition)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  select(-c("depth.factor","condition"))

# basic DESeq2 design
dds.gene = DESeqDataSetFromMatrix(countData = synthtest.counts,
                                  colData = synthtest.metadata,
                                  design = as.formula(~phenotype))
dds.gene.deg = DESeq(dds.gene)

#check data look reasonable
boxplot(log10(assays(dds.gene.deg)[["cooks"]]), range=0, las=2)
boxplot(log10(assays(dds.gene.deg)[["counts"]]), range=0, las=2)
plotDispEsts(dds.gene.deg)

# Okay that seems in the right ballpark for number of DEGs
table(results(dds.gene.deg)$padj<0.05)
thesedegs = row.names(subset(results(dds.gene.deg),padj<0.05))

# But are they the right genes? Here 1 = true DEG
truedegs = row.names(subset(synthtest@variable.annotations, differential.expression == 1))
table(thesedegs %in% truedegs) # Mostly - about 94%

# Check group means for a gene
mean(synthtest.counts[1,1:50])
mean(synthtest.counts[1,51:100])

# export so we can use for testing in sklearn
write.csv(synthtest.counts, file = "input/synthtest/synth_counts.csv")
write.csv(synthtest.metadata, file = "input/synthtest/synth_metadata.csv")
write.csv(synthtest@variable.annotations, file = "input/synthtest/synth_annotations.csv")


# Okay, we're actually having some issues over with sklearn 
# I suspect that the artificially low dispersions of non-DE genes relative to the strength of DE might be the issue?
# Try making all genes DE and then permuting the ones that we want to be non-DE

synthtest = generateSyntheticData(dataset = "testset_1", n.vars = 5000, 
                                  samples.per.cond = 50, n.diffexp = 5000, 
                                  repl.id = 1, seqdepth = 1e6, 
                                  fraction.upregulated = 0.5, 
                                  between.group.diffdisp = FALSE, 
                                  filter.threshold.total = 10, 
                                  filter.threshold.mediancpm = 0, 
                                  fraction.non.overdispersed = 0, 
                                  output.file = NULL)

# clean data
synthtest.counts = synthtest@count.matrix
synthtest.metadata = 
  synthtest@sample.annotations %>% 
  mutate(phenotype = as.factor(condition)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  select(-c("depth.factor","condition"))

# permute sample order for most of the genes
# permute function
permutegene = function(thisgene){
  permrow = unname(sample(synthtest.counts[thisgene,], replace = FALSE))
  return(permrow)
}
synthtest.permcounts = lapply(row.names(synthtest.counts), permutegene)
synthtest.permcounts = synthtest.permcounts %>% data.frame() %>% t() %>% 
  `row.names<-`(row.names(synthtest.counts)) %>%
  `colnames<-`(colnames(synthtest.counts))

#first 80% of rows become permuted
synthtest.counts.semiperm = synthtest.counts
synthtest.counts.semiperm[1:4900,] = synthtest.permcounts[1:4900,]

# check it worked 
# basic DESeq2 design
dds.gene = DESeqDataSetFromMatrix(countData = synthtest.counts.semiperm,
                                  colData = synthtest.metadata,
                                  design = as.formula(~phenotype))
dds.gene.deg = DESeq(dds.gene)

# A lot of false positives, it looks like? That maybe bodes well as a test set
table(results(dds.gene.deg)$padj<0.05)
thesedegs = row.names(subset(results(dds.gene.deg),padj<0.05))

# export so we can use for testing in sklearn
write.csv(synthtest.counts.semiperm, file = "input/synthtest/synth_counts_semiperm.csv")


# Read back in SVM weightings from sklearn
somecoefs = names(read.csv("testcoefs.csv", check.names = F))
testcoef = data.frame(name = row.names(synthtest.counts.semiperm), coef = as.numeric(somecoefs), num = 1:5000)
ggplot(testcoef, aes(x = num, y = coef)) +
  geom_point()

# Okay, so we do see that there are much stronger (more negative) coefs for the tail end of genes, which are the ones we left DE. But why are they all negative? That's weird, right? 
# Half the DE genes should be DE in each direction.... check that this is actually the case

