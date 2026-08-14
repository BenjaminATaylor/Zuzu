process ENET_BASIC {

  label 'mid_job'

  publishDir "$params.outdir/enet", pattern: 'enet_*'

  input:
  tuple path(samplesheet), path(countsframe)
  path library_sizes

  output:
  path "enet_gene_stability.csv", emit: table
  path "enet_stable_features.csv", emit: stable
  path "enet_summary.csv", emit: summary
  path "enet_fit_diagnostics.csv", emit: diagnostics
  path "enet_allocations_and_inner_folds.csv", emit: allocations
  path "enet_metadata.json", emit: metadata
  path "enet_top_features.png", emit: topFeaturesPlot
  path "enet_base_model_sizes.png", emit: baseModelSizesPlot
  path "enet_top_feature_selection_matrix.png", emit: selectionMatrixPlot

  script:
  task_seed = params.enet_seed.toInteger()
  """
  # Run the complete stability-selected elastic-net analysis
  zuzu_enet.R \
    --mode basic \
    --samplesheet "$samplesheet" \
    --counts "$countsframe" \
    --reflevel '${params.reflevel}' \
    --library-sizes "$library_sizes" \
    --alpha '${params.enet_alpha}' \
    --pairs '${params.enet_pairs}' \
    --pi '${params.enet_pi}' \
    --folds '${params.enet_folds}' \
    --nlambda '${params.enet_nlambda}' \
    --lambda-min-ratio '${params.enet_lambda_min_ratio}' \
    --seed '${task_seed}'
  """
}
