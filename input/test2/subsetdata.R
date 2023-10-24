setwd("~/Downloads")

samplesheet = read.csv(file = "bee_samplesheet_clean.csv")
samplesheet = subset(samplesheet, condition %in% c("Forager","Nurse")) %>%
  select(c("sample","condition")) %>%
  `colnames<-`(c("sample","phenotype"))

counts = read.delim("bee.salmon.merged.gene_counts.tsv") 

subcounts = counts %>%
  `colnames<-`(str_remove(colnames(counts),"X")) %>%
  select(c("gene_name",samplesheet$sample))

setwd("~/Repositories/Zuzu/input/test2/")
write.csv(samplesheet, row.names = FALSE, file = "samplesheet_phenotype.csv")
write.table(subcounts, sep = "\t", file = "salmon.merged.gene_counts.tsv")
