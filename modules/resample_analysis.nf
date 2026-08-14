process RESAMPLE_ANALYSIS {

  label 'mid_job'
  tag "${resampling_type}: target=${target}, repeat=${replicate}"

  input:
  tuple val(resampling_type), val(replicate), val(target), val(generation_seed)
  path samplesheet, stageAs: 'staged_inputs/source_samplesheet.csv'
  path countsframe, stageAs: 'staged_inputs/source_countsframe.tsv'
  val reflevel

  output:
  path "${resampling_type}_subsample_repeat_${replicate}_target_${target}_summary.csv", emit: summaries
  path "${resampling_type}_subsample_repeat_${replicate}_target_${target}_genes.csv", emit: genes
  path "${resampling_type}_subsample_repeat_${replicate}_target_${target}_allocations.csv", emit: allocations
  path "${resampling_type}_subsample_repeat_${replicate}_target_${target}_library_metrics.csv", emit: libraries

  script:
  prefix = "${resampling_type}_subsample_repeat_${replicate}_target_${target}"
  """
  # Run every enabled selector on one shared perturbed dataset
  zuzu_resample.R \
    --mode '${resampling_type}' \
    --repeat '${replicate}' \
    --target '${target}' \
    --seed '${generation_seed}' \
    --samplesheet '${samplesheet}' \
    --raw-counts '${countsframe}' \
    --reflevel '${reflevel}' \
    --output-prefix '${prefix}' \
    --enet '${params.enetstep}' \
    --enet-alpha '${params.enet_alpha}' \
    --enet-pairs '${params.enet_pairs}' \
    --enet-pi '${params.enet_pi}' \
    --enet-folds '${params.enet_folds}' \
    --enet-nlambda '${params.enet_nlambda}' \
    --enet-lambda-min-ratio '${params.enet_lambda_min_ratio}' \
    --enet-seed '${params.enet_seed}'
  """
}
