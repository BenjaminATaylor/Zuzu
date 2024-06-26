# Load a dataset for testing
library(tidyverse)
library(DESeq2)
library(compcodeR)
library(DescTools)

#### Load real data ####

setwd("~/Repositories/Zuzu/")
counts.true = read.csv("input/test5/salmon.merged.gene_counts.tsv", sep="\t", header=TRUE, row.names=1) %>%
  select(-c("gene_name")) %>% floor()
metadata.true = read.csv("input/test5/samplesheet.csv") %>% 
  mutate(phenotype = factor(phenotype, levels = c("nurse","forager")))

# check conformity
all(metadata.true$sample == colnames(counts.true))

# Extract basic information
ngenes = nrow(counts.true)
nsamples = ncol(counts.true)
ncounts = sum(counts.true)
samplesizes = table(metadata.true$phenotype)

#### Initial DESeq2 analysis ####

# Before moving onto simulation, we'll run DESeq2 for the real dataset
dds = DESeqDataSetFromMatrix(countData = counts.true,
                              colData = metadata.true,
                              design = ~ phenotype)
dds.true = DESeq(dds)
res.true = results(dds.true)
nDEGs = length(which(res.true$padj<0.05))
mcols(dds.true)


mcols(dds)$dispersion
# We'll also create individual DESeq2 objects using the group subsets, which will allow us to access e.g. per-group dispersions
dds.nurse = DESeqDataSetFromMatrix(countData = counts.true[,metadata.true$phenotype == "nurse"],
                                   colData = metadata.true %>% filter(phenotype == "nurse"),
                                   design = ~ 1)
dds.nurse = DESeq(dds.nurse)
dispersions.ref.true = mcols(dds.nurse)$dispersion 
dds.forager = DESeqDataSetFromMatrix(countData = counts.true[,metadata.true$phenotype == "forager"],
                                   colData = metadata.true %>% filter(phenotype == "forager"),
                                   design = ~ 1)
dds.forager = DESeq(dds.forager)
dispersions.alt.true = mcols(dds.forager)$dispersion 
# We do get some NAs here: we can see that these are all genes with uniform 0 counts in the given group
counts.true[which(is.na(mcols(dds.nurse)$dispersion )),metadata.true$phenotype == "nurse"]
counts.true[which(is.na(mcols(dds.forager)$dispersion )),metadata.true$phenotype == "forager"]
# Correct these to 0, since they do represent 0-dispersion genes in their respective groups
dispersions.ref.true[is.na(dispersions.ref.true)] = 0
dispersions.alt.true[is.na(dispersions.alt.true)] = 0


#### Synthetic data generation ####
# Generate a matching synthetic dataset using compcodeR
# Get means of the reference level in the real data
refmeans.true = rowMeans(counts.true[,metadata.true$phenotype == "nurse"])
dispersions.true = data.frame(ref = dispersions.ref.true, alt = dispersions.alt.true) %>% as.matrix()

# generate synthetic data
# TODO: adjust parameters based on the real input data 
data.compcoder = generateSyntheticData(dataset = "compcoder", repl.id = 123,
                                       n.vars = ngenes, 
                                       samples.per.cond = ceiling(nsamples/2),
                                       n.diffexp = nDEGs, 
                                       seqdepth = ncounts/nsamples, 
                                       relmeans = refmeans.true,
                                       dispersions = dispersions.true,
                                       fraction.upregulated = 0.5, 
                                       filter.threshold.total = 0, 
                                       filter.threshold.mediancpm = 0, 
                                       output.file = NULL)

counts.compcoder = data.frame(data.compcoder@count.matrix)
# Last thing to do here is to bring the sample size back in line, if compcodeR forced us to add samples to have even sample counts
keepsamples = row.names(rbind(subset(data.compcoder@sample.annotations,condition==1)[c(1:samplesizes[1]),],
                               subset(data.compcoder@sample.annotations,condition==2)[c(1:samplesizes[2]),]))
counts.compcoder = counts.compcoder[,keepsamples]
metadata.compcoder = select(data.compcoder@sample.annotations, condition) %>% 
  rownames_to_column() %>% `colnames<-`(c("sample","phenotype")) %>%
  subset(sample %in% keepsamples) %>%
  mutate(phenotype = fct_recode(factor(phenotype), "nurse" = "1", "forager" = "2")) 

# Run DESeq for the compcoder dataset
dds.compcoder = DESeqDataSetFromMatrix(countData = counts.compcoder,
                             colData = metadata.compcoder,
                             design = ~ phenotype)
dds.compcoder = DESeq(dds.compcoder)
res.compcoder = results(dds.compcoder)

#### Testing dataset conformity ####
#Now we can test how well this new dataset resembles the true dataset

#### First compare QQ plots of log-normalized values in the two datasets ####
quantiles = seq(0, 1, 0.00001)
qqdata.compcoder = data.frame(true = unlist(counts.true), compcoder = unlist(counts.compcoder))
ggplot(mapping = aes(x = quantile(log(qqdata.compcoder$true+1), quantiles, type = 5), 
                     y = quantile(log(qqdata.compcoder$compcoder+1), quantiles, type = 5))) + 
  geom_point() +
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2)  +
  geom_smooth(method = "lm") +
  coord_equal()

# Get goodness of fit metrics for the quantiles
ccc.compcoder.log = CCC(log(qqdata.compcoder$true+1), 
                        log(qqdata.compcoder$compcoder+1), 
                        ci = "z-transform", conf.level = 0.95)
ccc.compcoder.raw = CCC(qqdata.compcoder$true,
                        qqdata.compcoder$compcoder,
                        ci = "z-transform", conf.level = 0.95)
# Much higher for the log-transformed values
ccc.compcoder.log$rho.c
ccc.compcoder.raw$rho.c

# (Note that for these sample sizes, K-S tests are too sensitive to be of much use)

#### Next compare BCV vs dispersion plots ####
plotDispEsts(dds.true)
plotDispEsts(dds.compcoder)

# We can extract the parameters of the two dispersion functions, although I'm not sure that there's a statistical test we can apply here
dispersionFunction(dds.true)
dispersionFunction(dds.compcoder)


#### Next compare means vs variance plots ####
# mean of log2 expression and variance of log2 expression
rowVars = function(x) apply(x, 1, var)

meanvar.true = data.frame(mean = rowMeans(log2(counts.true+1)),
                         var = rowVars(log2(counts.true+1)))
meanvar.compcoder = data.frame(mean = rowMeans(log2(counts.compcoder+1)),
                               var = rowVars(log2(counts.compcoder+1)))

# Again, we can visually compare but extracting a statistical test is more difficult
ggplot(meanvar.true, aes(x = mean, y = var)) + 
  geom_point() +
  geom_smooth(method = "gam")
ggplot(meanvar.compcoder, aes(x = mean, y = var)) + 
  geom_point() +
  geom_smooth(method = "gam")

#### Next check if feature-feature correlations are preserved ####
# Take a subsample of genes, since doing this for all pairwise comparisons would be too computationally intense
varied.counts.true = counts.true[which(rowVars(log2(counts.true+1))!=0),] 
varied.counts.subset.true = varied.counts.true[sample(nrow(varied.counts.true), 25),]
# Get pairwise spearman cor between each pair of genes
corrs.mat.true = cor(t(varied.counts.subset.true), method = "spearman")
corrs.true = corrs.mat.true[lower.tri(corrs.mat.true)]

# Repeat for compcoder
varied.counts.compcoder = counts.compcoder[which(rowVars(log2(counts.compcoder+1))!=0),] 
varied.counts.subset.compcoder = varied.counts.compcoder[sample(nrow(varied.counts.compcoder), 25),]
# Get pairwise spearman cor between each pair of genes
corrs.mat.compcoder = cor(t(varied.counts.subset.compcoder), method = "spearman")
corrs.compcoder = corrs.mat.compcoder[lower.tri(corrs.mat.compcoder)]

corrs.ggdata = data.frame(true = corrs.true, compcoder = corrs.compcoder) %>% 
  reshape2::melt()
ggplot(corrs.ggdata, aes(x = variable, y = value)) + 
  geom_violin() +
  theme_bw()

# A simple t-test will do to compare these two distributions statistically
t.test(corrs.ggdata$value ~ corrs.ggdata$variable)


#### Volcano plots ####
# We can also compare the volcano plots of the two datasets

res.true %>% data.frame() %>%
  mutate(DEG = (padj<0.05 & abs(log2FoldChange)>1)) %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = DEG)) + 
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("black","red")) +
  geom_vline(xintercept = c(-1,1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  theme_bw()

res.compcoder %>% data.frame() %>%
  mutate(DEG = (padj<0.05 & abs(log2FoldChange)>1)) %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = DEG)) + 
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("black","red")) +
  geom_vline(xintercept = c(-1,1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  theme_bw()

#### Also generate PCAs for each method ####
counts_pca = function(counts, metadata){
  
  pca.counts = log2(counts+1)
  data.pca = prcomp(t(pca.counts))
  percent.var = (data.pca$sdev^2 / sum(data.pca$sdev^2))
  pca.out = list(values = data.frame(data.pca$x),
                 percent.var = percent.var)
  #connect to phenotypic data
  ggpcadata = pca.out$values %>%
    rownames_to_column(var = "sample") %>%
    left_join(metadata, by = "sample")
  
  print(pca.plot <- ggplot(ggpcadata, aes(x = PC1, y = PC2, color = phenotype)) +
          geom_point(size = 5, position = position_jitter(width = 0.5,height=0.5)) +
          #geom_text(vjust = -1) +
          xlab(paste0("PC",1,"( ",signif(pca.out$percent.var[1]*100, 3),"%)")) +
          ylab(paste0("PC",2,"( ",signif(pca.out$percent.var[2]*100, 3),"%)")) +
          theme_bw() +
          theme(aspect.ratio = 1,
                panel.grid = element_line(color = "grey95"),
                legend.text = element_text(size = 12),
                legend.title = element_text(face = "bold", size = 15),
                axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12),
                axis.title = element_text(face = "bold", size =12)))
}

counts_pca(counts.true, metadata.true)
counts_pca(counts.compcoder, metadata.compcoder)

#### Check per-gene expression level distributions ####
ggplot(data.frame(true = rowMeans(log(counts.true+1)), compcoder = rowMeans(log(counts.compcoder+1))) %>% 
         reshape2::melt(), aes(x = value, fill = variable)) + 
  geom_density(alpha = 0.3) +
  labs(x = "Per-gene log mean expression level", y = "Density", fill = "Method") +
  theme_bw()

#### Distribution of 0s per gene
ggplot(data.frame(true = rowSums(counts.true == 0), compcoder = rowSums(counts.compcoder == 0)) %>% 
         reshape2::melt(), aes(x = value, fill = variable)) + 
  geom_density(alpha = 0.3) +
  labs(x = "Per-gene 0 count frequency", y = "Density", fill = "Method") +
  theme_bw()


