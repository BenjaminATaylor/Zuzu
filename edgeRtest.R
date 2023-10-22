library(edgeR)

# Basic edgeR workflow, taken from the documentation
dds.edge = DGEList(counts=countsframe.clean,group=samplesheet$phenotype)
dds.edge = normLibSizes(dds.edge)
design.edge = model.matrix(~samplesheet$phenotype)
dds.edge.est = estimateDisp(dds.edge,design.edge)
edge.fit = glmQLFit(dds.edge.est,design.edge)
edge.qlf = glmQLFTest(edge.fit,coef=2)

# Tabular output (perform FDR correction directly rather than relying on topTags to do it for us)
edge.out = edge.qlf$table %>% mutate(FDR = p.adjust(PValue, method = "BH"))

write.csv(edge.out, row.names = T)


nDEGs = nrow(subset(edge.out,padj<0.05))
save(nDEGs,file = "nDEGs.RData")
save(dds.deg, file = "dds_deg.RData")