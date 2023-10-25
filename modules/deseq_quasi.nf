process DESEQ_QUASI {
  
  input: 
  tuple path(quasiframe), path(truedegs), path(samplesheet)

  output:
  path 'outframe.csv'

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse", quietly = TRUE)
  library("DESeq2", quietly = TRUE)
  library("pROC", quietly = TRUE)

  DEBUG=FALSE

  samplesheet = read.csv("$samplesheet")
  quasi.frame = read.csv("$quasiframe", row.names = 1)
  truedegs = get(load("$truedegs"))

  #generate model
  dds = DESeqDataSetFromMatrix(countData = quasi.frame,
                               colData = samplesheet,
                               design = as.formula(~phenotype))
  # Run the default analysis for DESeq2
  dds.deg.quasi = DESeq(dds, fitType = "parametric", betaPrior = FALSE)

  # 1. Power: the proportion of treu DEGs identified as DEGs in this comparison
  quasidegs = row.names(subset(results(dds.deg.quasi),padj<0.05))
  if(DEBUG){quasidegs = sample(row.names(results(dds.deg.quasi)),200)}
  #if(DEBUG){quasidegs = sample(truedegs, 70)}
  power = length(which(truedegs %in% quasidegs))/length(truedegs)

  # 2. FDP: the proportion of all DEGs that are false positives
  FDP = 1-(length(which(quasidegs %in% truedegs))/length(quasidegs))

  # 3. ROC AUC
  allgenes = row.names(results(dds.deg.quasi))
  truelabels = as.numeric(allgenes %in% truedegs)
  quasilabels = as.numeric(allgenes %in% quasidegs)
  AUC = auc(truelabels, quasilabels)
  # Remember that 0.5 = a truly random ROC, so for best effect we do (0.5-AUC)*2
  normAUC = (AUC-0.5)*2 # Greater deviation from 0 = better discrimination

  # Save outputs
  outframe = data.frame(power = power, FDP = FDP, normAUC = normAUC)
  write.csv(outframe, file="outframe.csv", row.names = FALSE)

  """

}