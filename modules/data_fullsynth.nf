process DATA_FULLSYNTH{

    label 'mid_job'
    tag "replicate=${x},ref=${refnum},alt=${altnum},depth=${depth}"

    input: 
    tuple val(x), val(refnum), val(altnum), val(depth)
    path(countsframe, stageAs: 'staged_inputs/source_countsframe.tsv')

    output:
    tuple val(depth),
            path("trueDEGs.RData"), 
            path("synthsheet.csv"), 
            path("synthcounts.csv"),
            val(refnum),
            val(altnum),
            val(x)

    script:
    """
    #!/usr/bin/env Rscript
    library(compcodeR)
    library(tidyverse)

    set.seed($x)
    refnum = as.integer($refnum)
    altnum = as.integer($altnum)
    depth = $depth
    simulation.size = max(refnum, altnum)
    stopifnot(refnum > 0L, altnum > 0L)

    # pull in original input data. We'll use this to create a synthetic dataset that's similar to the true dataset
    orig.countsframe = read.delim("$countsframe", row.names = 1, check.names = FALSE) %>% 
    select_if(is.numeric)
    
    # get number of genes, number of desired DEGs, and a vector of gene expression from real data
    ngenes = nrow(orig.countsframe)
    nsynthdegs = ceiling(ngenes*$params.synthdeg_pct)
    origmeans = rowMeans(orig.countsframe)
    
    # generate both conditions at the larger requested sample size
    # TODO: adjust parameters based on the real input data 
    synthdata = generateSyntheticData(dataset = "testset_1", n.vars = ngenes, 
                                    samples.per.cond = simulation.size, n.diffexp = nsynthdegs,
                                    repl.id = 1, seqdepth = depth, 
                                    relmeans = origmeans,
                                    fraction.upregulated = 0.5, 
                                    between.group.diffdisp = FALSE, 
                                    filter.threshold.total = 0,
                                    filter.threshold.mediancpm = 0, 
                                    fraction.non.overdispersed = 0, 
                                    output.file = NULL)

    # validate the generated sample annotations
    synthdata.metadata = synthdata@sample.annotations %>%
      rownames_to_column(var = "sample")
    generated.conditions = sort(unique(as.character(synthdata.metadata\$condition)))
    if(!identical(generated.conditions, c("1", "2"))){
      stop("compcodeR did not generate the expected condition codes 1 and 2")
    }
    if(anyDuplicated(synthdata.metadata\$sample) ||
       anyDuplicated(colnames(synthdata@count.matrix))){
      stop("Synthetic data contain duplicate sample identifiers")
    }
    if(!setequal(synthdata.metadata\$sample, colnames(synthdata@count.matrix))){
      stop("Synthetic count columns and sample annotations do not match")
    }

    # subset each condition to its exact requested size
    ref.samples = subset(synthdata.metadata, as.character(condition) == "1")\$sample
    alt.samples = subset(synthdata.metadata, as.character(condition) == "2")\$sample
    stopifnot(length(ref.samples) == simulation.size,
              length(alt.samples) == simulation.size)
    select.samples = function(sample.ids, requested.size){
      if(requested.size == length(sample.ids)){
        return(sample.ids)
      }
      sample.ids[sample.int(length(sample.ids), size = requested.size)]
    }
    chosen.samples = c(
      select.samples(ref.samples, refnum),
      select.samples(alt.samples, altnum)
    )
    synthdata.counts = synthdata@count.matrix[,chosen.samples,drop=FALSE]
    synthdata.metadata = synthdata.metadata[
      match(chosen.samples, synthdata.metadata\$sample),,
      drop=FALSE
    ]

    # apply the total-count filter after discarding excess samples
    keep.genes = rowSums(synthdata.counts) >= 10
    synthdata.counts = synthdata.counts[keep.genes,,drop=FALSE]
    synthdata.annotations = synthdata@variable.annotations[
      row.names(synthdata.counts),,
      drop=FALSE
    ]
    if(!identical(row.names(synthdata.counts), row.names(synthdata.annotations))){
      stop("Synthetic count rows and variable annotations do not match")
    }

    # relabel conditions and verify count-matrix alignment
    reflevel = "$params.reflevel"
    altlevel = if(reflevel == "nonref") "synthetic_alternative" else "nonref"
    synthdata.metadata\$phenotype = ifelse(
      as.character(synthdata.metadata\$condition) == "2",
      altlevel,
      reflevel
    )
    synthdata.metadata = select(synthdata.metadata, -c("depth.factor","condition"))
    stopifnot(nrow(synthdata.metadata) == refnum + altnum,
              ncol(synthdata.counts) == refnum + altnum,
              sum(synthdata.metadata\$phenotype == reflevel) == refnum,
              sum(synthdata.metadata\$phenotype == altlevel) == altnum,
              length(unique(synthdata.metadata\$phenotype)) == 2L,
              !anyDuplicated(chosen.samples),
              identical(colnames(synthdata.counts), synthdata.metadata\$sample))

    # extract names of true DEGs
    truedegs = row.names(subset(synthdata.annotations, differential.expression == 1))

    # Check that the synthetic truth universe is valid
    synthgenes = row.names(synthdata.counts)
    if(length(synthgenes) == 0){
      stop("Synthetic benchmark contains no genes after total-count filtering")
    }
    if(anyDuplicated(synthgenes)){
      stop("Synthetic benchmark contains duplicate gene identifiers")
    }
    if(anyDuplicated(truedegs)){
      stop("Synthetic truth set contains duplicate gene identifiers")
    }
    missing.truth = setdiff(truedegs, synthgenes)
    if(length(missing.truth) > 0){
      stop(sprintf(
        "Synthetic truth set contains %d genes absent from the generated count matrix",
        length(missing.truth)
      ))
    }
    if(length(truedegs) == 0){
      stop("Synthetic benchmark contains no true DEGs after total-count filtering")
    }
    if(length(truedegs) == length(synthgenes)){
      stop("Synthetic benchmark contains no background genes after total-count filtering")
    }

    # save outputs
    save(truedegs, file="trueDEGs.RData")
    write.csv(synthdata.metadata, file="synthsheet.csv", row.names = FALSE)
    write.csv(synthdata.counts, file="synthcounts.csv")
    """
}
