orig.countsframe = read.delim(file = "../input/test5/salmon.merged.gene_counts.tsv", row.names = 1, check.names = FALSE) %>% 
  select_if(is.numeric)

# What we want to do here is write a custom function that will take a dataframe of counts generate a synethetic dataset of the same size. Each gene in the synthetic dataset should follow the same distribution as the corresponding gene in the original dataset. However, in the sythetic dataset some genes should be differentially expressed between the two groups, while other should not be.

# The function should take the following arguments:
# - countsframe: a dataframe of counts where rows are genes and columns are samples
# - group: a vector of length equal to the number of columns in countsframe indicating the group of each sample
# - prop: the proportion of genes that should be differentially expressed
# - effect.size: the effect size of the differential expression
# - seed: a random seed

# The function should return a dataframe of counts where the rows are genes and the columns are samples

synthesise = function(countsframe, group, prop = 0.1, effect.size = 2, seed = 1){
  set.seed(seed)
  ngenes = nrow(countsframe)
  nsamples = ncol(countsframe)
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