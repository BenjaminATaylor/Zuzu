library(devtools)


ipak <- function(pkg, repository=c('CRAN', 'Bioconductor', 'github')){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  # new.pkg <- pkg
  if (length(new.pkg)) {
    if(repository=='CRAN') {
      install.packages(new.pkg, dependencies = TRUE)
    }
    if(repository=='Bioconductor') {
      if(strsplit(version[['version.string']], ' ')[[1]][3] > "3.6.0"){
        if (!requireNamespace("BiocManager")){
          install.packages("BiocManager")
        }
        BiocManager::install(new.pkg, dependencies=TRUE, ask=FALSE)
      }
      if(strsplit(version[['version.string']], ' ')[[1]][3] < "3.6.0"){
        stop(message("powsimR depends on packages that are only available in R 3.6.0 and higher."))
      }
    }
    if(repository=='github') {
      devtools::install_github(new.pkg, build_vignettes = FALSE, force = FALSE, dependencies=TRUE)
    }
  }
}

# CRAN PACKAGES
cranpackages <- c("broom", "cobs", "cowplot",
                  "data.table", "doParallel", "dplyr", "DrImpute",
                  "fastICA", "fitdistrplus", "foreach", "future",
                  "gamlss.dist", "ggplot2", "ggpubr", "ggstance", "grDevices",
                  "grid", "Hmisc", "kernlab", "MASS", "magrittr", "MBESS", "Matrix",
                  "matrixStats", "mclust", "methods", "minpack.lm", "moments", "msir",
                  "NBPSeq", "nonnest2", "parallel", "penalized", "plyr", "pscl",
                  "reshape2", "Rmagic", "rsvd", "Rtsne", "scales", "Seurat", "snow", "sctransform",
                  "stats", "tibble", "tidyr", "truncnorm", "VGAM", "ZIM", "zoo")
ipak(cranpackages, repository='CRAN')

# BIOCONDUCTOR
biocpackages <- c("bayNorm", "baySeq", "BiocGenerics", "BiocParallel",
                  "DESeq2", "EBSeq", "edgeR", "IHW", "iCOBRA",
                  "limma", "Linnorm", "MAST", "monocle", "NOISeq", "qvalue", "ROTS", "RUVSeq",
                  "S4Vectors", "scater", "scDD", "scde", "scone", "scran", "SCnorm",
                  "SingleCellExperiment", "SummarizedExperiment", "zinbwave")
ipak(biocpackages, repository='Bioconductor')

# GITHUB
githubpackages <- c('cz-ye/DECENT', 'nghiavtr/BPSC',
                    'mohuangx/SAVER', 'statOmics/zingeR',
                    'Vivianstats/scImpute', 'cran/Rmagic', 
                    'juba/rmdformats')
ipak(githubpackages, repository = 'github')



powsimRdeps <- data.frame(Package = c(cranpackages,
                                      biocpackages,
                                      sapply(strsplit(githubpackages, "/"), "[[", 2)),
                          stringsAsFactors = F)

ip <- as.data.frame(installed.packages()[,c(1,3:4)], stringsAsFactors = F)

ip.check <- cbind(powsimRdeps,
                  Version = ip[match(powsimRdeps$Package, rownames(ip)),"Version"])

table(is.na(ip.check$Version)) 

devtools::install_github("bvieth/powsimR", build_vignettes = TRUE, dependencies = FALSE, force = TRUE)
library(powsimR)











 data("CELseq2_Gene_UMI_Counts")
batch <- sapply(strsplit(colnames(CELseq2_Gene_UMI_Counts), "_"), "[[", 1)
Batches <- data.frame(Batch = batch,
                      stringsAsFactors = FALSE,
                      row.names = colnames(CELseq2_Gene_UMI_Counts))
data("GeneLengths_mm10")

# estimation
# powsimR can't handle underscores in the gene names properly, so we'll remove these temporarily then put them back in later
counts.true.tmp = data.frame(apply(counts.true, 2, as.integer))
rownames(counts.true.tmp) <- gsub("_", "---", rownames(counts.true))
# we only use samples from a single group for powsimR
counts.true.tmp = counts.true.tmp[,subset(metadata.true, phenotype == "nurse")$sample]
write.csv(counts.true.tmp,file = "~/Downloads/counts.true.tmp.csv")

estparam_gene <- estimateParam(countData = counts.true.tmp[c(1:50),c(10:19)],
                               readData = NULL,
                               batchData = NULL,
                               spikeData = NULL,
                               spikeInfo = NULL,
                               Lengths = NULL, 
                               MeanFragLengths = NULL,
                               RNAseq = 'bulk', Protocol = 'Read',
                               Distribution = 'NB', Normalisation = "MR",
                               GeneFilter = 0.1, SampleFilter = 3, 
                               sigma = 1.96, NCores = NULL, verbose = TRUE)

# Okay, there's some strange stuff going on here: this runs, but only with certain combinations of columns.
# I have raised an issue on github: https://github.com/bvieth/powsimR/issues/68






data("Bulk_Read_Counts")
data("GeneLengths_hg19")
estparam <- estimateParam(countData = Bulk_Read_Counts,
                          readData = NULL,
                          batchData = NULL,
                          spikeData = NULL,
                          spikeInfo = NULL,
                          Lengths = NULL,
                          MeanFragLengths = NULL,
                          RNAseq = 'bulk', Protocol = 'Read',
                          Distribution = 'NB', Normalisation = "MR",
                          GeneFilter = 0.1, SampleFilter = 3,
                          sigma = 1.96, NCores = NULL, verbose = TRUE)
