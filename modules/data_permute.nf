process DATA_PERMUTE{

  label 'small_job'
  tag "replicate=${x}"

  input: 
  tuple path(samplesheet), path(countsframe)
  val x

  output:
  tuple val(x), path("samplesheet_perm.csv"), path(countsframe), emit: frames
  tuple val(x), path("permutation_metadata.csv"), emit: metadata

  script:
  """
  #!/usr/bin/env Rscript
  library("tidyverse")

  samplesheet = read.csv("$samplesheet")

  # Permutation function (NB written with some help from colleagues)
  # In brief, how this works: we first find all possible fixed-size phenotype permutations, discard the ones that produce too large of a correlation with the original phenotype vector, and then sample uniformly from the remaining permutations
  permute_binary_labels = function(labels, reflevel, seed){
    labels = as.character(labels)
    stopifnot(!anyNA(labels))

    observed.levels = unique(labels)
    stopifnot(length(observed.levels) == 2L, reflevel %in% observed.levels)
    altlevel = setdiff(observed.levels, reflevel)
    stopifnot(length(altlevel) == 1L)

    original = as.integer(labels == altlevel)
    sample.count = length(original)
    alt.count = sum(original)
    ref.count = sample.count - alt.count
    stopifnot(alt.count > 0L, ref.count > 0L)

    # Enumerate attainable overlaps between the original and permuted alternative groups
    possible.overlaps = seq.int(max(0L, 2L * alt.count - sample.count), alt.count)
    numerator = as.double(sample.count) * as.double(possible.overlaps) -
                as.double(alt.count)^2
    denominator = as.double(alt.count) * as.double(ref.count)
    possible.correlations = numerator / denominator

    # For a threshold of e.g 0.6 = 3/5, use integer-equivalent arithmetic so as not to be affected by floating-point rounding.
    valid = 5 * abs(numerator) < 3 * denominator
    valid.overlaps = possible.overlaps[valid]
    if(length(valid.overlaps) == 0L){
      stop(sprintf(
        paste0(
          "No fixed-size phenotype permutation satisfies |correlation| < 0.6 ",
          "for group sizes %d and %d. Attainable correlations: %s"
        ),
        ref.count,
        alt.count,
        paste(signif(possible.correlations, 4), collapse = ", ")
      ))
    }

    # Weight overlaps by their number of possible assignments, making the draw
    # uniform over all label vectors that satisfy the correlation constraint.
    log.weights = lchoose(alt.count, valid.overlaps) +
                  lchoose(ref.count, alt.count - valid.overlaps)
    weights = exp(log.weights - max(log.weights))
    weights = weights / sum(weights)
    
    set.seed(as.integer(seed))
    selected.overlap = valid.overlaps[
      sample.int(length(valid.overlaps), 1L, prob = weights)
    ]

    take_positions = function(indices, size){
      if(size == 0L){
        integer()
      } else {
        indices[sample.int(length(indices), size, replace = FALSE)]
      }
    }

    new.alt.positions = c(
      take_positions(which(original == 1L), selected.overlap),
      take_positions(which(original == 0L), alt.count - selected.overlap)
    )
    permuted = integer(sample.count)
    permuted[new.alt.positions] = 1L
    permutation.correlation = cor(original, permuted)

    stopifnot(sum(permuted) == alt.count,
              is.finite(permutation.correlation),
              abs(permutation.correlation) < 0.6)

    list(
      labels = ifelse(permuted == 1L, altlevel, reflevel),
      correlation = permutation.correlation,
      reflevel = reflevel,
      altlevel = altlevel,
      ref.count = ref.count,
      alt.count = alt.count,
      alt.overlap = selected.overlap
    )
  }

  permutation = permute_binary_labels(
    samplesheet\$phenotype,
    reflevel = "$params.reflevel",
    seed = $x
  )

  samplesheet.perm = samplesheet
  samplesheet.perm\$phenotype = permutation\$labels
  non.phenotype.columns = setdiff(colnames(samplesheet), "phenotype")
  stopifnot(identical(samplesheet.perm[,non.phenotype.columns,drop=FALSE],
                      samplesheet[,non.phenotype.columns,drop=FALSE]))
  stopifnot(identical(table(samplesheet.perm\$phenotype),
                      table(as.character(samplesheet\$phenotype))))

  permutation.metadata = data.frame(
    replicate = as.integer($x),
    correlation = permutation\$correlation,
    max_abs_correlation = 0.6,
    reference = permutation\$reflevel,
    alternative = permutation\$altlevel,
    reference_count = permutation\$ref.count,
    alternative_count = permutation\$alt.count,
    alternative_overlap = permutation\$alt.overlap
  )

  write.csv(samplesheet.perm, "samplesheet_perm.csv", row.names = FALSE)
  write.csv(permutation.metadata, "permutation_metadata.csv", row.names = FALSE)

  """
}
