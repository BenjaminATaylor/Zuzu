# For this plotting, we need the following info for each type of analsis:

# How many DEGs were found in the original analysis? Just a number
# How many DEGs were found in each permutation (a list of numbers)

nDEGs = 100
deseq.permutes = get(load("/Users/benjamin/Repositories/Zuzu/work/7a/700cf3b3780c6e0a1f2e7c3786c6ce/nDEGs_list.RData"))

data.frame(analysis = "DESeq2", nDEGs = nDEGs, permutes = deseq.permutes)
data.frame(analysis = "DESeq2", nDEGs = nDEGs, permutes = deseq.permutes)


#Input data frame ends up looking like this:
permuteplot.input = data.frame(method = c(rep("DESeq",10),rep("edgeR",10)),
           nDEGs = c(rep(100,10),rep(130,10)),
           permutes = c(rnorm(10,100,5),rnorm(10,130,5)))

print(gg.permute <- ggplot(permuteplot.input, aes(x = method, y = permutes)) +
        geom_point(size = 3, alpha = 0.7) +
        geom_point(aes(x = method, y = nDEGs), 
                   size =4 , color = "black", fill = "red", shape = 23) +
        labs(x = "Method", y = "Number of DEGs") +
        theme_bw()
)

ggsave(gg.permute, 
       filename = "permute_plot.pdf",
       device = "pdf", bg = "transparent",
       width =  30, height = 20, units = "cm")



### Another plot we'd like to recreate is Fig 1B from Li et al. 
# This is a histogram using the fully permuted data, which shows the proportion of datasets in which a given gene was identified as DE when permuted (x axis) vs the number of genes for which that proportion was true (y axis)
# So to make this work, what we need is to get lists of all genes and whether they are DE or non-DE from a whole bunch of permutations. Example data:

foo = c("/Users/benjamin/Repositories/Zuzu/work/53/169bb981b7e290130fd97b23ad69c6/deseq_table.csv", "/Users/benjamin/Repositories/Zuzu/work/45/7efc46d023137303ded9715c7d117b/deseq_table.csv", "/Users/benjamin/Repositories/Zuzu/work/e8/9371b6fccf3fcfafa563dfe54ff1b9/deseq_table.csv")


thisfun = function(x){
  res <- read.csv(x,row.names = 1)
  df <- data.frame(result = (res$padj<0.05), row.names = row.names(res))
  return(df)
}

bar = lapply(foo,thisfun)
list_cbind(bar)

deseq.permdegs = lapply(foo,thisfun) %>% 
  list_cbind() %>%
  rowSums()

edger.permdegs = lapply(foo,thisfun) %>% 
  list_cbind() %>%
  rowSums() %>% sort()

permhist.in = cbind(deseq.permdegs,edger.permdegs) %>%
  `colnames<-`(c("DESeq2","edgeR")) %>%
  melt()

nperms = 100

# Generate breaks contingent on number of reps
cutbreaks = cut(seq(1,nperms+1),breaks = 5, right = TRUE) %>% levels() %>%
  str_match("\\([^,]*") %>%
  str_remove("\\(") %>%
  as.numeric() %>%
  round()
#deal with an idiosyncrasy that causes cut to return 0 as the start value of the first bin instead of 1
cutbreaks[1] = 1

if(DEBUG){permhist.in$value = blorg}

ggplot(permhist.in, aes(x = value)) +
  geom_histogram(color = "black", fill = "white", bins = 100) +
  scale_x_continuous(limits = c(1,(nreps+1)), breaks = cutbreaks) +
  labs(x = "Number of permutations in which gene was\nincorrectly identified as a DEG",
       y = "Number of genes") +
  theme_bw() +
  facet_grid(~Var2) +
  theme(strip.background = element_blank())








qux = data.frame(row.names = row.names(results(dds.deg)))
for(i in 1:10){
  
  qux = cbind(qux, data.frame(row.names = row.names(dds.deg), padj = runif(nrow(dds.deg),0,1)<0.25))  
  
}

nreps = 100


ggplot(data.frame(blorg), 
