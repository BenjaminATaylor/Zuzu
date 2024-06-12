orig.countsframe = read.delim(file = "../input/test5/salmon.merged.gene_counts.tsv", row.names = 1, check.names = FALSE) %>% 
  select_if(is.numeric)

orig.countsframe.clean = orig.countsframe %>% 
  filter(rowSums(.) > 10) %>% floor()

# for a given vector or numbers, we want to produce another vector of arbitrary length that draws from the same distribution as the input vector
# we want to do this in a way that generates new data, rather than simply sampling from the input vector
# we can do this by fitting a distribution to the input vector and then sampling from that distribution
library("MASS")

thisvector = rnbinom(n = 100, size = 1, mu = 10)
fitdistr(x = thisvector, densfun = "negative binomial")
ks.test(thisvector, "pnbinom", size = 1, mu = 10, simulate.p.value = T)

thisvector = unlist(orig.countsframe.clean[1,])


testnbinom = function(thisvector){
  thisdistr = fitdistr(x = thisvector,densfun = "negative binomial")
  ks.test(thisvector, "pnbinom", size = thisdistr$estimate["size"], 
          mu = thisdistr$estimate["mu"], simulate.p.value = T)
}

apply(orig.countsframe.clean[1:5,], 1, testnbinom)



# What we want to do here is write a custom function that will take a dataframe of counts generate a synethetic dataset of the same size. Each gene in the synthetic dataset should follow the same distribution as the corresponding gene in the original dataset. However, in the sythetic dataset some genes should be differentially expressed between the two groups, while other should not be.

# The function should take the following arguments:
# - countsframe: a dataframe of counts where rows are genes and columns are samples
# - group: a vector of length equal to the number of columns in countsframe indicating the group of each sample
# - prop: the proportion of genes that should be differentially expressed
# - effect.size: the effect size of the differential expression
# - seed: a random seed

# The function should return a dataframe of counts where the rows are genes and the columns are samples

synthesize = function(countsframe, group, prop = 0.1, effect.size = 2, seed = 1){
  set.seed(seed)
  ngenes = nrow(countsframe) ; nsamples = ncol(countsframe)
  ngenes.de = round(ngenes * prop)
  de.genes = sample(1:ngenes, ngenes.de)
  synthetic = countsframe
  for(i in 1:ngenes){
    if(i %in% de.genes){
      synthetic[i, group == 1] = countsframe[i, group == 1] + effect.size
    }
  }
  return(synthetic)
}