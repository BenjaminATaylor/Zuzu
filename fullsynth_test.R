library(compcodeR)
library(tidyverse)

setwd("/Users/benjamin/Repositories/Zuzu/")

set.seed(123)

refnum = 30
altnum = 30
depth = 1e7

# check which of the two sample sizes is larger
# we'll subset the columns for the other, since compcodeR cannot directly simulate uneven sample sizes
bignum = ifelse(altnum>refnum, altnum, refnum)
smallnum = ifelse(altnum>refnum, refnum, altnum)

bignum; smallnum

# generate synthetic data
# TODO: adjust parameters based on the real input data 
synthdata = generateSyntheticData(dataset = "qkwnqiuh", n.vars = 1000, 
                                  samples.per.cond = bignum, n.diffexp = 200, 
                                  repl.id = 1, seqdepth = depth, 
                                  fraction.upregulated = 0.5,
                                  #effect.size = c(rep(1,200),rep(2,800)),
                                  between.group.diffdisp = FALSE, 
                                  filter.threshold.total = 10, 
                                  filter.threshold.mediancpm = 0, 
                                  fraction.non.overdispersed = 0, 
                                  output.file = NULL)

params.reflevel = "foo"

# clean data, and change condition names to include our specified reference level
synthdata.counts = synthdata@count.matrix
synthdata.metadata = 
  synthdata@sample.annotations %>% 
  mutate(phenotype = as.factor(condition)) %>%
  mutate(phenotype = ifelse(phenotype == 2,"nonref",params.reflevel)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  select(-c("depth.factor","condition"))

# Plot to check that counts look the way we expect
synthdata.counts.scale = t(apply(synthdata.counts,1,scale)) %>% `colnames<-`(colnames(synthdata.counts)) 
absdiffs = rowMeans(synthdata.counts.scale[,1:30])-rowMeans(synthdata.counts.scale[,31:60])
ggplot(data.frame(diffs = absdiffs, index = 1:length(absdiffs)), aes(x = index, y = diffs)) +
  geom_point()

rowMeans(synthdata.counts.scale[,1:30])[500] - rowMeans(synthdata.counts.scale[,31:60])[500]


# extract names of true DEGs
truedegs = row.names(subset(synthdata@variable.annotations, differential.expression == 1))

# save outputs
save(depth, file="depth.RData")
save(truedegs, file="trueDEGs.RData")
write.csv(synthdata.metadata, file="input/synthtest/synth_sheet.csv")
write.csv(synthdata.counts, file="input/synthtest/synth_counts.csv")

# basic DESeq2 design
dds.gene = DESeqDataSetFromMatrix(countData = synthdata.counts,
                                  colData = synthdata.metadata,
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
truedegs = row.names(subset(synthdata@variable.annotations, differential.expression == 1))
table(thesedegs %in% truedegs) # Mostly - about 94%
table(truedegs %in% thesedegs)

# Check group means for a gene
mean(synthdata.counts[1,1:30])
mean(synthdata.counts[1,31:60])

# export so we can use for testing in sklearn
write.csv(synthdata.counts, file = "input/synthtest/synth_counts.csv")
write.csv(synthdata.metadata, file = "input/synthtest/synth_metadata.csv")
write.csv(synthdata@variable.annotations, file = "input/synthtest/synth_annotations.csv")

# Okay, we're actually having some issues over with sklearn 
# I suspect that the artificially low dispersions of non-DE genes relative to the strength of DE might be the issue?
# Try making all genes DE and then permuting the ones that we want to be non-DE

synthdata = generateSyntheticData(dataset = "testset_1", n.vars = 5000, 
                                  samples.per.cond = 50, n.diffexp = 5000, 
                                  repl.id = 1, seqdepth = 1e6, 
                                  fraction.upregulated = 0.5, 
                                  between.group.diffdisp = FALSE, 
                                  filter.threshold.total = 10, 
                                  filter.threshold.mediancpm = 0, 
                                  fraction.non.overdispersed = 0, 
                                  output.file = NULL)

# clean data
synthdata.counts = synthdata@count.matrix
synthdata.metadata = 
  synthdata@sample.annotations %>% 
  mutate(phenotype = as.factor(condition)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  select(-c("depth.factor","condition"))

# permute sample order for most of the genes
# permute function
permutegene = function(thisgene){
  permrow = unname(sample(synthdata.counts[thisgene,], replace = FALSE))
  return(permrow)
}
synthdata.permcounts = lapply(row.names(synthdata.counts), permutegene)
synthdata.permcounts = synthdata.permcounts %>% data.frame() %>% t() %>% 
  `row.names<-`(row.names(synthdata.counts)) %>%
  `colnames<-`(colnames(synthdata.counts))

#first 80% of rows become permuted
synthdata.counts.semiperm = synthdata.counts
synthdata.counts.semiperm[1:4900,] = synthdata.permcounts[1:4900,]

# check it worked 
# basic DESeq2 design
dds.gene = DESeqDataSetFromMatrix(countData = synthdata.counts.semiperm,
                                  colData = synthdata.metadata,
                                  design = as.formula(~phenotype))
dds.gene.deg = DESeq(dds.gene)

# A lot of false positives, it looks like? That maybe bodes well as a test set
table(results(dds.gene.deg)$padj<0.05)
thesedegs = row.names(subset(results(dds.gene.deg),padj<0.05))

# export so we can use for testing in sklearn
write.csv(synthdata.counts.semiperm, file = "input/synthtest/synth_counts_semiperm.csv")


# Read back in SVM weightings from sklearn
somecoefs = names(read.csv("testcoefs.csv", check.names = F))
testcoef = data.frame(name = row.names(synthdata.counts.semiperm), coef = as.numeric(somecoefs), num = 1:5000)
ggplot(testcoef, aes(x = num, y = coef)) +
  geom_point()

# Okay, so we do see that there are much stronger (more negative) coefs for the tail end of genes, which are the ones we left DE. But why are they all negative? That's weird, right? 
# Half the DE genes should be DE in each direction.... check that this is actually the case
synthdata.counts.DE = synthdata.counts[4901:5000,]
table(rowMeans(synthdata.counts.DE[1:50,]) > rowMeans(synthdata.counts.DE[51:100,]))
# Yep!



# Try again for synth data: do values and weightings match to what we expect?
somecoefs = names(read.csv("testcoefs_new.csv", check.names = F))
somemetadata = read.csv("testmetdata_new.csv", check.names = F, header = F)
somevalues = read.csv("testdata_new.csv", check.names = F, header = F) %>%
  t() %>% `colnames<-`(somemetadata$V1) %>% data.frame() %>% rownames_to_column(var = "gene")

plotframe = reshape2::melt(somevalues) %>% mutate(variable = substr(variable,2,2))

# scatter plot is not very informative
ggplot(plotframe, aes(x = variable, y = value)) +
  geom_point() +
  facet_grid(~gene)
# make_classification begins by putting all informative feature sin the first set of variables, so we should have the first 10 here informative and the remaining random
ggplot(plotframe, aes(x = variable, y = value)) +
  geom_boxplot() +
  facet_grid(~gene, scales = "free")

# Hmm, with that amount of noise it's not surprising we didn't find much.

# Let's try again with a larger class separation
# Try again for synth data: do values and weightings match to what we expect?
somecoefs = names(read.csv("testcoefs_new_double.csv", check.names = F))
somemetadata = read.csv("testmetdata_new_double.csv", check.names = F, header = F)
somevalues = read.csv("testdata_new_double.csv", check.names = F, header = F) %>%
  t() %>% `colnames<-`(somemetadata$V1) %>% data.frame() %>% rownames_to_column(var = "gene")

plotframe = reshape2::melt(somevalues) %>% mutate(variable = substr(variable,2,2))

ggplot(plotframe, aes(x = variable, y = value)) +
  geom_boxplot() +
  facet_grid(~gene, scales = "free")

# Ooooookay so comparing this to the previous plots, we can immediately see at least two issues:
# 1. There are only 3 informative features, not the expected 10
# 2. Those three features are not the first three, but rather seem to be 4, 12 and 22, scattered among the others (this latter part shouldn't have affected classification, but is further indication that make_classification is working in an unexpected manner.)

# are these at least the features that are being picked out by the classifier?
plot(somecoefs)
which(abs(as.numeric(somecoefs))>1)
# Yep, those are the right ones at least


# Some more testing with random make_classification setups
somemetadata = read.csv("testmetdata_new_thistest.csv", check.names = F, header = F)
somevalues = read.csv("testdata_new_thistest.csv", check.names = F, header = F) %>%
  t() %>% `colnames<-`(somemetadata$V1) %>% data.frame() %>% rownames_to_column(var = "gene")

plotframe = reshape2::melt(somevalues) %>% mutate(variable = substr(variable,2,2))

ggplot(plotframe, aes(x = variable, y = value)) +
  geom_boxplot() +
  facet_grid(~gene, scales = "free")



# index0 = which(substr(colnames(somevalues),2,2)=="0")                      
# index1 = which(substr(colnames(somevalues),2,2)=="1")     
# 
# somevals0 = somevalues[,index0]
# somevals1 = somevalues[,index1]





