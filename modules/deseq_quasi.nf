process DESEQ_QUASI {
  
  label 'mid_job'
  tag "replicate=${x},ref=${refnum},alt=${altnum}"
  
  input: 
  tuple path(quasiframe),
        path(truedegs),
        path(samplesheet),
        val(x),
        val(refnum),
        val(altnum)

  output:
  path 'outframe.csv'

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse", quietly = TRUE)
  library("DESeq2", quietly = TRUE)

  DEBUG=FALSE

  samplesheet = read.csv("$samplesheet")
  quasi.frame = read.csv("$quasiframe", row.names = 1, check.names = FALSE)
  truedegs = get(load("$truedegs"))
  replicate.id = $x
  ref.num = as.integer($refnum)
  alt.num = as.integer($altnum)
  reflevel = "$params.reflevel"

  # The quasi sample sheet has crossed another CSV boundary; restore its levels.
  samplesheet\$phenotype = relevel(factor(samplesheet\$phenotype), ref = reflevel)
  stopifnot(nlevels(samplesheet\$phenotype) == 2L,
            identical(levels(samplesheet\$phenotype)[1], reflevel))
  altlevel = setdiff(levels(samplesheet\$phenotype), reflevel)
  stopifnot(length(altlevel) == 1L)

  observed.refnum = sum(samplesheet\$phenotype == reflevel)
  observed.altnum = sum(samplesheet\$phenotype != reflevel)
  stopifnot(observed.refnum == ref.num)
  stopifnot(observed.altnum == alt.num)
  stopifnot(nrow(samplesheet) == ref.num + alt.num)
  stopifnot(ncol(quasi.frame) == ref.num + alt.num)
  stopifnot(identical(colnames(quasi.frame), samplesheet\$sample))

  # Generate output metrics only if there are at least 10 'true' DEGs, otherwise outputs NAs
  if(length(truedegs)>10){
    #generate model
    dds = DESeqDataSetFromMatrix(countData = quasi.frame,
                                colData = samplesheet,
                                design = as.formula(~phenotype))
    # Run the default analysis for DESeq2
    dds.deg.quasi = DESeq(dds, fitType = "parametric", betaPrior = FALSE)
    deseq.results = results(dds.deg.quasi,
                            contrast = c("phenotype", altlevel, reflevel))

    # 1. Power: the proportion of treu DEGs identified as DEGs in this comparison
    quasidegs = row.names(subset(deseq.results,padj<0.05))
    if(DEBUG){quasidegs = sample(row.names(deseq.results),200)}
    #if(DEBUG){quasidegs = sample(truedegs, 70)}
    power = length(which(truedegs %in% quasidegs))/length(truedegs)

    # 2. FDP: the proportion of all DEGs that are false positives
    FDP = 1-(length(which(quasidegs %in% truedegs))/length(quasidegs))

    # 3. Youden's index at the DEG threshold
    allgenes = row.names(deseq.results)
    truth = allgenes %in% truedegs
    called = allgenes %in% quasidegs
    sensitivity = sum(called & truth) / sum(truth)
    false_positive_rate = sum(called & !truth) / sum(!truth)
    youden_index = sensitivity - false_positive_rate

  } else {
    power = NA
    FDP = NA
    youden_index = NA
  }

  outframe = data.frame(power = power,
                        FDP = FDP,
                        youden_index = youden_index,
                        replicate = replicate.id,
                        refnum = ref.num,
                        altnum = alt.num,
                        totalnum = nrow(samplesheet))
  write.csv(outframe, file="outframe.csv", row.names = FALSE)
  """

}
