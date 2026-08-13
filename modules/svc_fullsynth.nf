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
  path "svc_fullsynth_*_n_supported_features.txt", emit: nSupported
  path "svc_fullsynth_*_gate.csv", emit: gate
  path "svc_fullsynth_*_gate_null.csv", emit: gateNull
  path "svc_fullsynth_*_gate.png", emit: gatePlot

  script:
  depth_token = depth.toLong()
  output_prefix = "svc_fullsynth_depth_${depth_token}_ref_${refnum}_alt_${altnum}_replicate_${x}"
  """
  # Fit adaptive SVC RFE-CV without changing the existing full-synthetic dataset
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
    --coarse-grid '${params.svc_coarse_grid}' \
    --refinement-points '${params.svc_refinement_points}' \
    --n-jobs '${params.njobs}' \
    --gate true \
    --gate-permutations '${params.svc_gate_permutations}' \
    --gate-alpha '${params.svc_gate_alpha}' \
    --gate-seed '${params.svc_gate_seed}' \
    --gate-early-stop '${params.svc_gate_early_stop}' \
    --analysis-type fullsynthetic \
    --replicate '${x}' \
    --depth '${depth}' \
    --reference-samples '${refnum}' \
    --alternative-samples '${altnum}' \
    --output-prefix '${output_prefix}'
  """
}
