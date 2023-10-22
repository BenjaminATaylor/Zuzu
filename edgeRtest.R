library(edgeR)

# Basic edgeR workflow, taken from the documentation
dds.edge = DGEList(counts=countsframe.clean,group=samplesheet$phenotype)
dds.edge = normLibSizes(dds.edge)
design.edge = model.matrix(~samplesheet$phenotype)
dds.edge.est = estimateDisp(dds.edge,design.edge)
edge.fit = glmQLFit(dds.edge.est,design.edge)
edge.qlf = glmQLFTest(edge.fit,coef=2)

# Tabular output (rather than topTags, which only gives the most highly-DE results)
edge.out = edge.qlf$table
write.csv(edge.out, row.names = T)
