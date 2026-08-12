process CLEANINPUTS {

    label 'small_job'
    debug true

    input: 
    path(samplesheet, stageAs: 'staged_inputs/source_samplesheet.csv')
    path(countsframe, stageAs: 'staged_inputs/source_countsframe.tsv')
    val reflevel

    output:
    tuple path("samplesheet.csv"), path("countsframe_clean.csv")

    script:
    """
    #!/usr/bin/env Rscript
    library("tidyverse", quietly = TRUE)

    # Setup
    samplesheet = read.csv("$samplesheet")
    countsframe = read.delim("$countsframe", row.names = 1, check.names = FALSE) %>% 
        select_if(is.numeric)

    # Check count and phenotype matrix conformity
    stopifnot(identical(sort(samplesheet\$sample), sort(colnames(countsframe))))
    # Order both frames
    samplesheet = samplesheet[order(samplesheet\$sample),]
    countsframe = countsframe[,order(colnames(countsframe))]

    # Check that the 'phenotype' conditions has exactly two conditions, and set a reference level
    stopifnot("$reflevel" %in% levels(as.factor(samplesheet\$phenotype)))
    samplesheet\$phenotype = as.factor(samplesheet\$phenotype) %>% relevel(ref = "$reflevel")
    stopifnot(length(levels(samplesheet\$phenotype))==2)

    # Round all values to integers
    countsframe = round(countsframe)

    # Retain genes with at least 1 count per million in at least the smaller
    # group size of samples; group membership itself does not affect which
    # genes pass.
    count.matrix = as.matrix(countsframe)
    if(any(!is.finite(count.matrix))){
      stop("Count matrix contains missing or non-finite values")
    }
    if(any(count.matrix < 0)){
      stop("Count matrix contains negative values")
    }
    library.sizes = colSums(count.matrix)
    if(any(library.sizes <= 0)){
      stop("CPM filtering requires every sample to have a positive library size")
    }

    minimum.samples = min(as.integer(table(samplesheet\$phenotype)))
    at.least.one.cpm = sweep(
      count.matrix,
      MARGIN = 2L,
      STATS = library.sizes / 1e6,
      FUN = ">="
    )
    keep.rows = rowSums(at.least.one.cpm) >= minimum.samples
    filtered.genes = sum(!keep.rows)
    filtered.proportion = filtered.genes / length(keep.rows)
    if(filtered.proportion > 0.5){
      filter.warning = sprintf(
        paste0(
          "WARNING [CLEANINPUTS]: the common CPM/prevalence filter removed ",
          "%d of %d genes (%.1f%%), exceeding the 50%% warning threshold. ",
          "Check library sizes and whether the 1 CPM in %d samples rule is ",
          "appropriate for this dataset."
        ),
        filtered.genes,
        length(keep.rows),
        100 * filtered.proportion,
        minimum.samples
      )
      cat(filter.warning, "\n")
      warning(filter.warning, call. = FALSE)
    }
    if(!any(keep.rows)){
      stop("No genes pass the common CPM/prevalence filter")
    }
    countsframe.clean = countsframe[keep.rows,,drop=FALSE]

    write.csv(samplesheet,"samplesheet.csv", row.names = FALSE)
    write.csv(countsframe.clean,"countsframe_clean.csv")
    """
}
