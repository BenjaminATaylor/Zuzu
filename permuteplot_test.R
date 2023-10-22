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
                   size =3, color = "red") +
        labs(x = "Method", y = "Number of DEGs") +
        theme_bw()
)

