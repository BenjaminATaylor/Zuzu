library("seqgendiff")
library(airway)
library(tidyverse)
data(airway)
colData(airway)

library(limma)
library(sva)
library(qvalue)

# Treat the true signal in the data as a surrogate variable, which we'll account for having added in a new variable as the one to test
true_sv <- model.matrix(~phenotype, data = metadata.true)[, -1]

# Here we apply a signal to 10% of genes (1-prop_null), and with a log effect distributed in a normal distribution with mean 0 and standard deviation 0.8
# The ditribution of effect sizes could be tweaked to fit more closely to the true data, perhaps by trying to fit a distribution to the effect sizes of DEGs in the true data
thout <- thin_2group(mat = as.matrix(counts.true), 
                     prop_null = 0.9, 
                     signal_fun = stats::rnorm,
                     signal_params = list(mean = 0, sd = 0.8))




X <- cbind(thout$design_obs, thout$designmat)
Y <- log2(thout$mat + 0.5)
n_sv <- num.sv(dat = Y, mod = X)
svout <- sva(dat = Y, mod = X, n.sv = n_sv)

vout <- voom(counts = thout$mat, design = cbind(X, svout$sv))
lout <- lmFit(vout)
eout <- eBayes(lout)
qout <- qvalue(p = eout$p.value[, 2])
bhat <- eout$coefficients[, 2]

plot(thout$coefmat, 
     bhat, 
     xlab = "True Coefficients", 
     ylab = "Estimated Coefficients")
abline(0, 1, col = 2, lwd = 2)

is_null_gene <- abs(thout$coefmat) < 10^-6
boxplot(qout$qvalues ~ is_null_gene,
        xlab = "Null Gene",
        ylab = "q-value")













#OLD
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
