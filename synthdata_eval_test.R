# Load a dataset for testing
library(tidyverse)
library(DESeq2)
library(compcodeR)
library(DescTools)
library(SimSeq)

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
## Generate a matching synthetic dataset using compcodeR
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

## Generate another synthetic dataset using seqgendiff
# Here we apply a signal to 10% of genes (1-prop_null), and with a log effect distributed in a normal distribution with mean 0 and standard deviation 0.8
prop.deg = nDEGs/nrow(counts.true)
counts.ref.true = counts.true[,metadata.true$phenotype == "nurse"]

# The distribution of effect sizes could be tweaked to fit more closely to the true data, perhaps by trying to fit a distribution to the effect sizes of DEGs in the true data
thout = thin_2group(mat = as.matrix(counts.ref.true),
                    prop_null = 1-prop.deg, 
                    signal_fun = stats::rnorm,
                    signal_params = list(mean = 0, sd = 0.8))
thout$mat

# Hmm but here we have 2 issues. 1) We only have the sample size from the control data, not the original data;
# and 2) The effects have been distributed to a random subset of samples, instead of us choosing which groups are true and false
# 2) is fine (can rearrange posthoc) but how to deal with 1)?
# For now we'll just have to live with it
#Extract counts
counts.seqgendiff = thout$mat %>% 
  data.frame() %>%
  `colnames<-`(paste0("sample_",1:ncol(.))) %>%
  `rownames<-`(row.names(counts.true.nurse))
# Extract metadata
tmp = thout$designmat %>% 
  data.frame() %>% 
  mutate(phenotype = factor(ifelse(P1 == 0, "nurse", "forager"), levels = c("nurse","forager"))) %>%
  mutate(sample = paste0("sample_",1:nrow(.))) %>%
  select(sample, phenotype)
metadata.seqgendiff = tmp[order(tmp$phenotype),]
counts.seqgendiff = counts.seqgendiff[,metadata.seqgendiff$sample]
all(colnames(counts.seqgendiff) == metadata.seqgendiff$sample)

# Run DESeq for the seqgendiff dataset
dds.seqgendiff = DESeqDataSetFromMatrix(countData = counts.seqgendiff,
                                       colData = metadata.seqgendiff,
                                       design = ~ phenotype)
dds.seqgendiff = DESeq(dds.seqgendiff)
res.seqgendiff = results(dds.seqgendiff)

## Generate another synthetic dataset using simSeq
# Note an issue here- due to the way the simulation is set up, we can only simulate half the sample size of the real data using simSeq. For now this is just something we'll have to live with
data.simseq = SimData(counts = counts.true, 
                      treatment = metadata.true$phenotype,
                      sort.method = "unpaired",
                      switch.trt = FALSE,
                      k.ind = floor(min(table(metadata.true$phenotype))/2),  # = samples per treatment group
                      n.genes = nrow(counts.true), 
                      n.diff = floor(nrow(counts.true)*0.1))
#Get metadata
metadata.simseq = data.simseq$treatment %>% 
  data.frame() %>% 
  mutate(phenotype = ifelse(. == 0, "nurse", "forager")) %>%
  mutate(sample = paste0("sample_",row_number())) %>%
  select(sample, phenotype)
# Get counts data
counts.simseq = data.simseq$counts %>% `colnames<-`(metadata.simseq$sample)
all(colnames(counts.simseq)==metadata.simseq$sample)

# Run DESeq for the simseq dataset
dds.simseq = DESeqDataSetFromMatrix(countData = counts.simseq,
                                        colData = metadata.simseq,
                                        design = ~ phenotype)
dds.simseq = DESeq(dds.simseq)
res.simseq = results(dds.simseq)



#### Testing dataset conformity ####
#Now we can test how well these new datasets resemble the true dataset
methods = c("compcoder","seqgendiff","simseq")

#### First compare QQ plots of log-normalized values in the datasets ####
quantiles = seq(0, 1, 0.00001)

for(method in methods){
  # Get the objects for the method
  dds = get(paste0("dds.",method))
  counts = get(paste0("counts.",method))
  metadata = get(paste0("metadata.",method))
  res = get(paste0("res.",method))
  
  qq.true = quantile(log(unlist(counts.true)+1), quantiles, type = 5)
  qq.sim = quantile(log(unlist(counts)+1), quantiles, type = 5)
  
  ggplot(mapping = aes(x = qq.true, 
                       y = qq.sim)) + 
    geom_point() +
    geom_abline(aes(slope = 1, intercept = 0), linetype = 2)  +
    geom_smooth(method = "lm") +
    coord_equal() +
    theme_bw()
  
  # Get goodness of fit metrics for the quantiles
  ccc.log = CCC(log(qq.true+1), log(qq.sim+1), ci = "z-transform", conf.level = 0.95)
  ccc.raw = CCC(qq.true, qq.sim, ci = "z-transform", conf.level = 0.95)
  ccc.log$rho.c
  ccc.raw$rho.c
  
  out = data.frame(method = method, ccclog = ccc.log$rho.c$est, cccraw = ccc.raw$rho.c$est)
  
  if(method == methods[1]){
    qqout = out
  } else {
    qqout = rbind(qqout, out)
  }
  
  #pseudocode: we can either plot the qqs together,or make single plot with the qq fits for each method plotted    next to each other (ideally we would have many fits from many repeated simulations, so then we could plot as a   box or something like that)
}
# (Note that for these sample sizes, K-S tests are too sensitive to be of much use)

#### Next compare BCV vs dispersion plots ####
plotDispEsts2 <- function(dds){
  as.data.frame(mcols(dds)) %>% 
    select(baseMean, dispGeneEst, dispFit, dispersion) %>%
    melt(id.vars="baseMean") %>% 
    filter(baseMean>0) %>% 
    ggplot(aes(x=baseMean, y=value, colour=variable)) + 
    geom_point(size=0.1) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_bw() + 
    ylab("Dispersion") + 
    xlab("BaseMean") +
    scale_colour_manual(
      values=c("Black", "#e41a1c", "#377eb8"), 
      breaks=c("dispGeneEst", "dispFit", "dispersion"), 
      labels=c("Estimate", "Fit", "Final"),
      name=""
    ) +
    guides(colour = guide_legend(override.aes = list(size=2)))
}

mcols.true = as.data.frame(mcols(dds)) %>% 
  select(baseMean, dispGeneEst, dispFit, dispersion) %>%
  mutate(method = "true")

for(method in methods){
  
  dds = get(paste0("dds.",method))
  mcols = as.data.frame(mcols(dds)) %>% 
    select(baseMean, dispGeneEst, dispFit, dispersion) %>%
    mutate(method = method)
  
  if(method == methods[1]){
    mcolsout = rbind(mcols.true, mcols)
  } else {
    mcolsout = rbind(mcolsout, mcols)
  }
  
}

mcolsout%>% 
  melt(id.vars=c("baseMean","method")) %>% 
  mutate(method = factor(method, levels = c("true",methods))) %>%
  filter(baseMean>0) %>% 
  ggplot(aes(x=baseMean, y=value, colour=variable)) + 
  geom_point(size=0.1) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_bw() + 
  ylab("Dispersion") + 
  xlab("BaseMean") +
  scale_colour_manual(
    values=c("Black", "#e41a1c", "#377eb8"), 
    breaks=c("dispGeneEst", "dispFit", "dispersion"), 
    labels=c("Estimate", "Fit", "Final"),
    name=""
  ) +
  facet_wrap(~method) +
  guides(colour = guide_legend(override.aes = list(size=2)))

# We can extract the parameters of the two dispersion functions, although I'm not sure that there's a statistical test we can apply here
dispersionFunction(dds.true)
dispersionFunction(dds.compcoder)

#### Next compare means vs variance plots ####
# mean of log2 expression and variance of log2 expression
rowVars = function(x) apply(x, 1, var)

meanvar.true = data.frame(mean = rowMeans(log2(counts.true+1)),
                         var = rowVars(log2(counts.true+1)),
                         method = "true")

for(method in methods){
  
  counts = get(paste0("counts.",method))
  meanvar = data.frame(mean = rowMeans(log2(counts+1)),
                       var = rowVars(log2(counts+1)),
                       method = method)
  
  if(method == methods[1]){
    meanvarout = rbind(meanvar.true, meanvar)
  } else {
    meanvarout = rbind(meanvarout, meanvar)
  }
  
}

meanvarout = mutate(meanvarout, method = factor(method, levels = c("true",methods)))

# Again, we can visually compare but extracting a statistical test is more difficult
ggplot(meanvarout, aes(x = mean, y = var)) +
  geom_point() +
  geom_smooth(method = "gam") +
  facet_wrap(~method)


#### Next check if feature-feature correlations are preserved ####
# Take a subsample of genes, since doing this for all pairwise comparisons would be too computationally intense
varied.counts.true = counts.true[which(rowVars(log2(counts.true+1))!=0),] 
varied.counts.subset.true = varied.counts.true[sample(nrow(varied.counts.true), 500),]
# Get pairwise spearman cor between each pair of genes
corrs.mat.true = cor(t(varied.counts.subset.true), method = "spearman")
corrs.true = corrs.mat.true[lower.tri(corrs.mat.true)]

# Generalise across methods
for(method in methods){
  
  # Get the objects for the method
  counts = get(paste0("counts.",method))
  
  # Repeat for compcoder
  varied.counts = counts[which(rowVars(log2(counts+1))!=0),] 
  varied.counts.subset = varied.counts[sample(nrow(varied.counts), 500),]
  # Get pairwise spearman cor between each pair of genes
  corrs.mat = cor(t(varied.counts.subset), method = "spearman")
  corrs = corrs.mat[lower.tri(corrs.mat)]
  
  corrs.ggdata = data.frame(true = corrs.true, compcoder = corrs) %>% 
    reshape2::melt()
  ggplot(corrs.ggdata, aes(x = variable, y = value)) + 
    geom_violin() +
    theme_bw()

  # A simple t-test and also K-S test
  t.test(corrs.ggdata$value ~ corrs.ggdata$variable)
  ks.test(corrs.true, corrs)
  # A better was of testing this might be the Wasserstein distance?
  library("transport")
  wasserstein1d(a = corrs.true, b = corrs)
  
  out = data.frame(method = method, wasserstein = wasserstein1d(a = corrs.true, b = corrs))
  
  if(method == methods[1]){
    corrsout = out
  } else {
    corrsout = rbind(corrsout, out)
  }
}

corrsout



#### Volcano plots ####
# We can also compare the volcano plots of the datasets

volcdat.true <- res.true %>% data.frame() %>%
  mutate(DEG = (padj<0.05 & abs(log2FoldChange)>1)) %>% 
  mutate(method = "true")

for(method in methods){
  
  volcdat <- get(paste0("res.",method)) %>% data.frame() %>%
    mutate(DEG = (padj<0.05 & abs(log2FoldChange)>1)) %>%
    mutate(method = method)
  
  if(method == methods[1]){
    volcdats = rbind(volcdat.true, volcdat)
  } else {
    volcdats = rbind(volcdats, volcdat)
  }
}
#plot
volcdats = mutate(volcdats, method = factor(method, levels = c("true",methods)))
ggplot(volcdats,aes(x = log2FoldChange, y = -log10(padj), color = DEG)) + 
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("black","red")) +
  geom_vline(xintercept = c(-1,1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  theme_bw() +
  facet_wrap(~method)

#### Also generate PCAs for each method ####
data_pca = function(counts, metadata, method = "true"){
  
  pca.counts = log2(counts+1)
  data.pca = prcomp(t(pca.counts))
  percent.var = (data.pca$sdev^2 / sum(data.pca$sdev^2))
  pca.out = list(values = data.frame(data.pca$x),
                 percent.var = percent.var)
  #connect to phenotypic data
  ggpcadata = pca.out$values %>%
    select(paste0("PC",c(1:4))) %>%
    rownames_to_column(var = "sample") %>%
    left_join(metadata, by = "sample") %>%
    mutate(method = method)
  
  return(ggpcadata)
}

for(method in c("true",methods)){
  
  ggpcadata = data_pca(get(paste0("counts.",method)), get(paste0("metadata.",method)), method)
  
  if(method == "true"){
    ggpcadataout = ggpcadata
  } else {
    ggpcadataout = rbind(ggpcadataout, ggpcadata)
  }
  
}

# reorder for plotting
ggpcadataout = mutate(ggpcadataout, method = factor(method, levels = c("true",methods)))
# Plot
print(pca.plot <- ggplot(ggpcadataout, aes(x = PC1, y = PC2, color = phenotype)) +
        geom_point(size = 5, position = position_jitter(width = 0.5,height=0.5)) +
        #geom_text(vjust = -1) +
        xlab("PC1") + ylab("PC2") +
        theme_bw() +
        facet_wrap(~method) +
        theme(aspect.ratio = 1,
              panel.grid = element_line(color = "grey95"),
              legend.text = element_text(size = 12),
              legend.title = element_text(face = "bold", size = 15),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              axis.title = element_text(face = "bold", size =12)))

#### Check per-gene expression level distributions ####
logs.true = data.frame(value = rowMeans(log(counts.true+1)), method = "true")
for(method in methods){
  counts = get(paste0("counts.",method))
  logs = data.frame(value = rowMeans(log(counts+1)), method = method)
  
  # A simple t-test and also K-S test
  t.test(logs.true$value, logs$value)
  suppressWarnings(ks.test(logs.true$value, logs$value))
  out = data.frame(method = method, wasserstein = wasserstein1d(a = logs.true$value, b = logs$value))

  if(method == methods[1]){
    logsout = rbind(logs, logs.true)
    logcompsout = out
  } else {
    logsout = rbind(logsout, logs)
    logcompsout = rbind(logcompsout, out)
  }
}

logsout = mutate(logsout, method = factor(method, levels = c("true",methods)))
ggplot(logsout, aes(x = value, fill = "grey90")) +
  geom_density(alpha = 0.3) +
  scale_fill_discrete(guide = "none") +
  labs(x = "Per-gene log mean expression level", y = "Density") +
  theme_bw() +
  facet_wrap(~method)
logcompsout


#### Distribution of 0s per gene
zeros.true = data.frame(value = rowSums(counts.true == 0), method = "true")

for(method in methods){
  counts = get(paste0("counts.",method))
  zeros = data.frame(value = rowSums(counts == 0), method = method)

  out = data.frame(method = method, wasserstein = wasserstein1d(a = zeros.true$value, b = zeros$value))
  
  if(method == methods[1]){
    zerosout = rbind(zeros, zeros.true)
    zeroscompsout = out
  } else {
    zerosout = rbind(zerosout, zeros)
    zeroscompsout = rbind(zeroscompsout, out)
  }
}

zerosout = mutate(zerosout, method = factor(method, levels = c("true",methods)))
ggplot(zerosout, aes(x = value, fill = "grey90")) + 
  geom_density(alpha = 0.3) +
  scale_fill_discrete(guide = "none") +
  labs(x = "Per-gene 0 count frequency", y = "Density", fill = "Method") +
  theme_bw() +
  facet_wrap(~method)
zeroscompsout



#### Distribution of per-gene variances
vars.true = data.frame(value = log(meanvar.true$var+1), method = "true")

for(method in methods){
  
  counts = get(paste0("counts.",method))
  meanvar = data.frame(mean = rowMeans(log2(counts+1)),
                       var = rowVars(log2(counts+1)))
  vars = data.frame(value = log(meanvar$var+1), method = method)
  
  out = data.frame(method = method, wasserstein = wasserstein1d(a = vars.true$value, b = vars$value))
  
  if(method == methods[1]){
    varsout = rbind(vars, vars.true)
    varscompsout = out
  } else {
    varsout = rbind(varsout, vars)
    varscompsout = rbind(varscompsout, out)
  }

}

varsout = mutate(varsout, method = factor(method, levels = c("true",methods)))
ggplot(varsout, aes(x = value, fill = "grey50")) +
  scale_fill_discrete(guide = "none") +
  geom_density(alpha = 0.3) +
  labs(x = "Per-gene log variance", y = "Density", fill = "Method") +
  theme_bw() +
  facet_wrap(~method)

