process QUALITYCONTROL {

  publishDir "$params.outdir"

  input: 
  tuple path(samplesheet), path(countsframe)

  output:
  path "QC_plots.pdf"

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse")

  # Here we're looking to generate some nice diagnostic output plots
  # Some of this draws from this standard protocol: https://f1000research.com/articles/4-1070

  library("DESeq2")
  library("pheatmap")
  library("RColorBrewer")
  library("reshape2")
  library("ggpubr")
  library("ggplotify")

  samplesheet = read.csv("$samplesheet")
  countsframe.clean = read.csv("$countsframe",row.names = 1)

  vst.counts = countsframe.clean %>% as.matrix() %>% varianceStabilizingTransformation()

  ## Heatmap of sample distances
  sampleDists = dist(t(vst.counts))
  sampleDistMatrix = as.matrix(sampleDists)
  #rownames(sampleDistMatrix) <- paste( rld\$dex, rld\$cell, sep="-" )
  #colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  print(heatmap <- pheatmap(sampleDistMatrix,
                            clustering_distance_rows=sampleDists,
                            clustering_distance_cols=sampleDists,
                            annotation_col = column_to_rownames(samplesheet,var = "sample"),
                            annotation_row = column_to_rownames(samplesheet,var = "sample"),
                            col=colors))

  ## PCA, with annotations
  #create pca object
  data.pca = prcomp(t(vst.counts))
  #extract PC data
  percent.var = (data.pca\$sdev^2 / sum(data.pca\$sdev^2))
  pca.out = list(values = data.frame(data.pca\$x),percent.var = percent.var)
  #connect to phenotypic data
  ggpcadata = pca.out\$values %>%
    rownames_to_column(var = "sample") %>%
    left_join(samplesheet,
              by = "sample")
  #plot
  print(gg.pca <- ggplot(ggpcadata, aes(x = PC1, y = PC2, color = phenotype, label = sample)) +
          geom_point(size = 5) +
          xlab(paste0("PC",1,"( ",signif(pca.out\$percent.var[1]*100, 3),"%)")) +
          ylab(paste0("PC",2,"( ",signif(pca.out\$percent.var[2]*100, 3),"%)")) +
          theme_bw() +
          theme(aspect.ratio = 1,
                text = element_text(size = 12),
                title = element_text(face = "bold", size = 15)))

  ## DESeq2 diagnostic plots are good here too:
  #generate a DESeq2 model
  dds.gene = DESeqDataSetFromMatrix(countData = countsframe.clean,
                                    colData = samplesheet,
                                    design = as.formula(~ phenotype))
  dds.gene.deg.local = DESeq(dds.gene, fitType = "local", betaPrior = FALSE)
  plotDispEsts(dds.gene.deg.local)
  dds.gene.deg.parametric = DESeq(dds.gene, fitType = "parametric", betaPrior = FALSE)
  plotDispEsts(dds.gene.deg.parametric)
  # Check for outliers with cook's distance
  print(gg.cooks <- ggplot(melt(log10(assays(dds.gene.deg.parametric)[["cooks"]]), na.rm = T), aes(x=Var2, y=value)) +
          geom_boxplot() +
          labs(x = "Sample", y = "Cook's distance") +
          theme_bw() +
          theme(axis.text.x = element_text(angle=90)))
  # And with counts,including pseudocount to accommodate logging
  print(gg.counts <- ggplot(melt(log10(assays(dds.gene.deg.parametric)[["counts"]]+1), na.rm = T), aes(x=Var2, y=value)) +
          geom_boxplot() +
          labs(x = "Sample", y = "Counts") +
          theme_bw() +
          theme(axis.text.x = element_text(angle=90)))
  ggout = ggarrange(as.ggplot(heatmap), gg.pca, gg.cooks, gg.counts)
  ggsave(ggout, 
         filename = "QC_plots.pdf",
         device = "pdf", bg = "transparent",
         width =  40, height = 40, units = "cm")

  """

}