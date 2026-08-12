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

  script:
  """
  # Fit the baseline SVC-RFECV procedure using pre-filter library sizes
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
    --step '${params.stepsize}' \
    --n-jobs '${params.njobs}' \
    --analysis-type basic \
    --output-prefix svc_basic
  """
}
