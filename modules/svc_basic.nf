process SVC_BASIC {

  label 'big_job'

  publishDir "$params.outdir/svctables", pattern: 'svc_basic_*'

  input:
  tuple path(samplesheet), path(countsframe)
  path library_sizes
  path library_genes

  output:
  path "svc_basic_table.csv", emit: table
  path "svc_basic_curve.csv", emit: curve
  path "svc_basic_folds.csv", emit: folds
  path "svc_basic_fs.png", emit: plot
  path "svc_basic_model.pickle", emit: model
  path "svc_basic_metadata.json", emit: metadata
  path "svc_basic_n_selected_features.txt", emit: nFeatures
  path "svc_basic_n_supported_features.txt", emit: nSupported
  path "svc_basic_gate.csv", emit: gate
  path "svc_basic_gate_null.csv", emit: gateNull
  path "svc_basic_gate.png", emit: gatePlot

  script:
  """
  # Fit the baseline adaptive SVC RFE-CV procedure using pre-filter library sizes
  zuzu_svc.py \
    --samplesheet "$samplesheet" \
    --counts "$countsframe" \
    --library-sizes "$library_sizes" \
    --library-genes "$library_genes" \
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
    --analysis-type basic \
    --output-prefix svc_basic
  """
}
