tcounts = read.delim("/Users/benjamin/Downloads/TCGA-LIHC.htseq_counts.tsv", row.names = 1)
tcounts = (tcounts[!str_detect(row.names(tcounts),pattern = "__"),])

tdata = read.delim("/Users/benjamin/Downloads/TCGA-LIHC.GDC_phenotype.tsv") %>% 
  select(c("submitter_id.samples","sample_type.samples")) %>%
  `colnames<-`(c("sample","phenotype")) %>%
  mutate(sample = str_replace_all(sample,"-",".")) %>%
  subset(sample %in% colnames(tcounts))

tdata$phenotype = ifelse(str_detect(tdata$phenotype,"Normal"),"Normal","Tumor")

dim(tcounts[,subset(tdata, phenotype == "Normal")$sample])
dim(tcounts[,subset(tdata, phenotype == "Tumor")$sample])

groupsize = 30
picksamples = c(subset(tdata, phenotype == "Tumor")$sample[1:groupsize],
                subset(tdata, phenotype == "Normal")$sample[1:groupsize])

genesize = nrow(tcounts)
pickgenes = rownames(tcounts)[1:genesize]
testcounts = tcounts[pickgenes,picksamples]
tesdata = tdata[match(picksamples,tdata$sample),]
table(testdata$sample == colnames(testcounts))

testcounts = floor(2^(testcounts-1))
#cut genes with few counts
testcounts = testcounts[(rowSums(testcounts)>=ncol(testcounts)),]

write.csv(testcounts, file = "input/test4/tumorcleancounts.csv")
write.csv(tesdata, row.names = F, file = "input/test4/tumorcleandata.csv")


# basic DESeq2 design
dds.gene = DESeqDataSetFromMatrix(countData = testcounts,
                                  colData = testdata,
                                  design = as.formula(~phenotype))
dds.gene.deg = DESeq(dds.gene)

#check data look reasonable
boxplot(log10(assays(dds.gene.deg)[["cooks"]]), range=0, las=2)
boxplot(log10(assays(dds.gene.deg)[["counts"]]), range=0, las=2)
plotDispEsts(dds.gene.deg)

# Okay that seems in the right ballpark for number of DEGs
table(results(dds.gene.deg)$padj<0.05)
thesedegs = row.names(subset(results(dds.gene.deg),padj<0.05))


     