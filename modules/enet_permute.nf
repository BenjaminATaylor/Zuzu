process ENET_PERMUTE {

  label 'mid_job'
  tag "replicate=${x}"

  publishDir "$params.outdir/enet", pattern: 'enet_permute_*'

  input:
  tuple val(x), path(samplesheet), path(countsframe)
  path library_sizes

  output:
  tuple val(x),
        path("enet_permute_replicate_${x}_stable_features.csv"),
        emit: stable
  tuple val(x),
        path("enet_permute_replicate_${x}_summary.csv"),
        emit: summary
  tuple val(x),
        path("enet_permute_replicate_${x}_fit_diagnostics.csv"),
        path("enet_permute_replicate_${x}_metadata.json"),
        emit: diagnostics

  script:
  task_seed = params.enet_seed.toInteger() + x.toInteger()
  """
  # Run the complete selector on one shared phenotype-label permutation
  zuzu_enet.R \
    --mode compact \
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

  # Give every retained outer-permutation artifact a unique name
  mv enet_stable_features.csv \
    "enet_permute_replicate_${x}_stable_features.csv"
  mv enet_summary.csv \
    "enet_permute_replicate_${x}_summary.csv"
  mv enet_fit_diagnostics.csv \
    "enet_permute_replicate_${x}_fit_diagnostics.csv"
  mv enet_metadata.json \
    "enet_permute_replicate_${x}_metadata.json"
  """
}
