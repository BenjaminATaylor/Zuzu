library(SPsimSeq) 

# load the Zhang bulk RNA-seq data (availabl with the package) 
data("zhang.data.sub") 
# filter genes with sufficient expression (important step to avoid bugs) 
zhang.counts <- zhang.data.sub$counts # Despite the implication of the documentation, these are indeed counts
MYCN.status  <- zhang.data.sub$MYCN.status #The grouping variable

set.seed(123) #Set seed for reproducibility
# simulate data
sim.data.bulk <- SPsimSeq(n.sim = 1, s.data = zhang.counts,
                          group = MYCN.status, n.genes = 3000, batch.config = 1,
                          group.config = c(0.5, 0.5), tot.samples = ncol(zhang.counts), 
                          pDE = 0.1, lfc.thrld = 0.5, result.format = "list", return.details = TRUE)

sim.data.bulk1 <- sim.data.bulk$sim.data.list[[1]]
head(sim.data.bulk1$counts[, seq_len(5)])  # count data

head(sim.data.bulk1$colData)        # sample info
head(sim.data.bulk1$rowData)        # gene info

geneDens = evaluateDensities(sim.data.bulk, newData = rownames(zhang.counts)[1])
#This returns for every sample, the midpoints (mids) and associated densities (gy)

# compare the distributions of the mean expressions, variability, 
# and fraction of zero counts per gene
library(LSD) # for generating heatmap plots
# normalize counts for comparison  
Y0.log.cpm <- log2(edgeR::cpm(zhang.counts)+1)
Y1.log.cpm <- log2(edgeR::cpm(sim.data.bulk1$counts)+1)
Y0.log.cpm <- Y0.log.cpm[rowMeans(Y0.log.cpm>0)>=0.1, ]
Y1.log.cpm <- Y1.log.cpm[rowMeans(Y1.log.cpm>0)>=0.1, ]
rowVars <- function(X){apply(X, 1, var, na.rm=TRUE)}
rowCVs <- function(X){apply(X, 1, function(x) sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE))}
par(mfrow=c(1, 3))
boxplot(list(real.data=log(colSums(zhang.counts)), 
             simulated.data=log(sim.data.bulk1$colData$sim.Lib.Size)), 
        main="library size") 
boxplot(list(real.data=rowMeans(Y0.log.cpm), 
             simulated.data=rowMeans(Y1.log.cpm)), 
        main="mean expression of genes") 
boxplot(list(real.data=rowVars(Y0.log.cpm), 
             simulated.data=rowVars(Y1.log.cpm)), 
        main="variance of gene expressions") 

# compare the relationship between the mean and variability
par(mfrow=c(1,3), mar=c(4,4,4,1))
heatscatter(rowMeans(Y0.log.cpm), rowCVs(Y0.log.cpm), ylim=c(0, 6), xlim=c(0, 16),
            colpal="bl2gr2rd", main="real data", xlab="mean log2-CPM", 
            ylab="coefficients of variation", cexplot=0.5, alpha = 60, cex.lab=1.25)
heatscatter(rowMeans(Y1.log.cpm), rowCVs(Y1.log.cpm), ylim=c(0, 6), xlim=c(0, 16),
            main="SPsimSeq", xlab="mean log2-CPM", ylab="coefficients of variation", 
            cexplot=0.5, alpha = 60, colpal="bl2gr2rd", cex.lab=1.25)

n.gride <- 1000
min.g   <- seq(0, 20, length.out = n.gride+1)[-n.gride]
max.g   <- seq(0, 20, length.out = n.gride+1)[-1] 
mid.g   <- (min.g+max.g)/2
f.real  <- vapply(seq_len(n.gride), FUN.VALUE = double(1), function(r){
  x <- Y0.log.cpm[rowMeans(Y0.log.cpm)<=max.g[r] & rowMeans(Y0.log.cpm)>min.g[r],]
  y <- ifelse(!is.null(dim(x)), mean(rowCVs(x)), mean(sd(x)/mean(x))) 
  y
})
f.SPsim <- vapply(seq_len(n.gride), FUN.VALUE = double(1), function(r){
  x <- Y1.log.cpm[rowMeans(Y1.log.cpm)<=max.g[r] & rowMeans(Y1.log.cpm)>min.g[r],]
  y <- ifelse(!is.null(dim(x)), mean(rowCVs(x)), mean(sd(x)/mean(x))) 
  y
})

sm1 <- loess(I(f.SPsim-f.real)~mid.g) 
plot(mid.g, f.SPsim-f.real, xlim=c(0, 14), col="lightskyblue", pch=20, cex.lab=1.25,
     cex.main=1.4, main="SPsimSeq - real data", ylab="difference", xlab="mean log2-CPM")
lines(mid.g,predict(sm1, newdata = mid.g), col="blue", lwd=3) 


# compare the correlation between genes and samples 
cor.mat.Y0 <- cor(t(Y0.log.cpm))
cor.mat.Y1 <- cor(t(Y1.log.cpm)) 
cor.vec.Y0 <- cor.mat.Y0[upper.tri(cor.mat.Y0)]
cor.vec.Y1 <- cor.mat.Y1[upper.tri(cor.mat.Y1)] 
par(mfrow=c(1,3), mar=c(4,4,3.5,1))
hist(cor.vec.Y0, nclass = 30, probability = TRUE, 
     border="gray", col="steelblue1", main="real data", xlab="Genewise correlations", 
     ylim=c(0, 3.5), xlim=c(-1, 1), cex.lab=1.25)
hist(cor.vec.Y1, nclass = 30, probability = TRUE, border="gray",
     col="steelblue1",  main="SPsimSeq", xlab="Genewise correlations",
     ylim=c(0, 3.5), xlim=c(-1, 1), cex.lab=1.25)
plot(seq(-1, 1, 0.1), seq(-1, 1, 0.1), type="n", xlab="quantile (real data)", 
     ylab="quantile (simulated data)",  main="correlation quantile-quantile plot")
abline(0, 1, col="gray")
points(quantile(cor.vec.Y0, seq(0, 1, 0.001)), quantile(cor.vec.Y1, seq(0, 1, 0.001)), 
       col="blue", pch=20, cex=1.5, cex.lab=1.25)         



#### For our data ####
# To make SPsimSeq work, we need to ensure even sample sizes and remove genes with zero counts
evensample = ncol(counts.tmp) - (ncol(counts.tmp) %% 2)
counts.tmp = counts.true[-c(which(rowSums(counts.true) == 0)),]
spsimseq = SPsimSeq(n.sim = 1, 
                    s.data = as.matrix(counts.tmp),
                    group = metadata.true$phenotype, 
                    n.genes = nrow(counts.tmp),
                    batch.config = 1,
                    group.config = c(0.5, 0.5), 
                    tot.samples = evensample, 
                    pDE = 0.1, 
                    lfc.thrld = 0.5, 
                    result.format = "list", 
                    return.details = TRUE)

