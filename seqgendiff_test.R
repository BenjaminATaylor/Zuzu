library("seqgendiff")
library(airway)
library(tidyverse)
data(airway)
colData(airway)

# take just control data as input
counts.true.nurse = counts.true[, metadata.true$phenotype == "nurse"]

# Here we apply a signal to 10% of genes (1-prop_null), and with a log effect distributed in a normal distribution with mean 0 and standard deviation 0.8
# The ditribution of effect sizes could be tweaked to fit more closely to the true data, perhaps by trying to fit a distribution to the effect sizes of DEGs in the true data
thout = thin_2group(mat = as.matrix(counts.true.nurse), 
            prop_null = 0.9, 
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
  mutate(phenotype = ifelse(P1 == 0, "nurse", "forager")) %>%
  mutate(sample = paste0("sample_",1:nrow(.))) %>%
  select(sample, phenotype)
metadata.seqgendiff = tmp[order(tmp$phenotype),]
counts.seqgendiff = counts.seqgendiff[,metadata.seqgendiff$sample]
all(colnames(counts.seqgendiff) == metadata.seqgendiff$sample)
