process SVC_FULLSYNTH {

  label 'big_job'
  tag "replicate=${x},ref=${refnum},alt=${altnum},depth=${depth}"

  publishDir "$params.outdir/svctables", pattern: 'svc_fullsynth_*'

  input:
  tuple val(depth),
        path(truedegs),
        path(samplesheet),
        path(countsframe),
        val(refnum),
        val(altnum),
        val(x)

  output:
  tuple val(depth),
        path(truedegs),
        val(refnum),
        val(altnum),
        val(x),
        path("svc_fullsynth_*_table.csv"), emit: result
  path "svc_fullsynth_*_curve.csv", emit: curve
  path "svc_fullsynth_*_folds.csv", emit: folds
  path "svc_fullsynth_*_fs.png", emit: plot
  path "svc_fullsynth_*_model.pickle", emit: model
  path "svc_fullsynth_*_metadata.json", emit: metadata
  path "svc_fullsynth_*_n_selected_features.txt", emit: nFeatures

  script:
  depth_token = depth.toLong()
  output_prefix = "svc_fullsynth_depth_${depth_token}_ref_${refnum}_alt_${altnum}_replicate_${x}"
  """
  # Fit SVC-RFECV without changing the existing full-synthetic dataset
  zuzu_svc.py \
    --samplesheet "$samplesheet" \
    --counts "$countsframe" \
    --reference '${params.reflevel}' \
    --pare '${params.pareml}' \
    --svc-c '${params.svc_c}' \
    --class-weight '${params.svc_class_weight}' \
    --folds '${params.svc_folds}' \
    --scoring '${params.svc_scoring}' \
    --seed '${params.svc_seed}' \
    --step '${params.stepsize}' \
    --n-jobs '${params.njobs}' \
    --analysis-type fullsynthetic \
    --replicate '${x}' \
    --depth '${depth}' \
    --reference-samples '${refnum}' \
    --alternative-samples '${altnum}' \
    --output-prefix '${output_prefix}'
  """
}
