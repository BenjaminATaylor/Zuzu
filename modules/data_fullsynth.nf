process DATA_FULLSYNTH{

    input: 
    tuple val(x), val(refnum), val(altnum), val(depth)

    output:
    tuple val(depth),
            path("trueDEGs.RData"), 
            path("synthsheet.csv"), 
            path("synthcounts.csv"),
            val(refnum),
            val(x)

    script:
    """
    #!/usr/bin/env Rscript
    library(compcodeR)
    library(tidyverse)

    set.seed($x)
    refnum = $refnum
    depth = $depth

    # if we wanted to, we could include uneven sample sizes- this to be added in future if desired
    # check which of the two sample sizes is larger
    # we'll subset the columns for the other, since compcodeR cannot directly simulate uneven sample sizes
    #bignum = ifelse(altnum>refnum, altnum, refnum)
    #smallnum = ifelse(altnum>refnum, refnum, altnum)

    # generate synthetic data
    # TODO: adjust parameters based on the real input data 
    synthdata = generateSyntheticData(dataset = "testset_1", n.vars = 2000, 
                                    samples.per.cond = refnum, n.diffexp = 200, 
                                    repl.id = 1, seqdepth = depth, 
                                    fraction.upregulated = 0.5, 
                                    between.group.diffdisp = FALSE, 
                                    filter.threshold.total = 10, 
                                    filter.threshold.mediancpm = 0, 
                                    fraction.non.overdispersed = 0, 
                                    output.file = NULL)

    # Clean data, and change condition names to include our specified reference level
    synthdata.counts = synthdata@count.matrix
    synthdata.metadata = synthdata@sample.annotations %>% 
    mutate(phenotype = as.factor(condition)) %>%
    mutate(phenotype = ifelse(phenotype == 2,"nonref","$params.reflevel")) %>%
    rownames_to_column(var = "sample") %>%
    select(-c("depth.factor","condition"))

    # extract names of true DEGs
    truedegs = row.names(subset(synthdata@variable.annotations, differential.expression == 1))

    # save outputs
    save(truedegs, file="trueDEGs.RData")
    write.csv(synthdata.metadata, file="synthsheet.csv")
    write.csv(synthdata.counts, file="synthcounts.csv")
    """
}
