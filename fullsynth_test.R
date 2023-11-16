library(compcodeR)
library(dplyr)
library(DESeq2)

setwd("Repositories/Zuzu/")

# generate synthetic data with 10000 genes, 20 samples/group, and 1200 differentially expressed genes
synthtest = generateSyntheticData(dataset = "testset_1", n.vars = 10000, 
                                  samples.per.cond = 20, n.diffexp = 1200, 
                                  repl.id = 1, seqdepth = 1e7, 
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
  mutate(condition = as.factor(condition)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  select(-c("depth.factor"))

# basic DESeq2 design
dds.gene = DESeqDataSetFromMatrix(countData = synthtest.counts,
                                  colData = synthtest.metadata,
                                  design = as.formula(~condition))
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

# export so we can use for testing in sklearn
write.csv(synthtest.counts, file = "input/synthtest/synth_counts.csv")
write.csv(synthtest.metadata, file = "input/synthtest/synth_metadata.csv")
write.csv(synthtest@variable.annotations, file = "input/synthtest/synth_annotations.csv")
