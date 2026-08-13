process SVC_PERMUTE {

  label 'big_job'
  tag "replicate=${x}"

  publishDir "$params.outdir/svctables", pattern: 'svc_permute_*'

  input:
  tuple val(x), path(samplesheet), path(countsframe)
  path library_sizes
  path library_genes

  output:
  path "svc_permute_replicate_${x}_table.csv", emit: outfile
  path "svc_permute_replicate_${x}_curve.csv", emit: curve
  path "svc_permute_replicate_${x}_folds.csv", emit: folds
  path "svc_permute_replicate_${x}_fs.png", emit: plot
  path "svc_permute_replicate_${x}_model.pickle", emit: model
  path "svc_permute_replicate_${x}_metadata.json", emit: metadata
  path "svc_permute_replicate_${x}_n_selected_features.txt", emit: nFeatures
  path "svc_permute_replicate_${x}_n_supported_features.txt", emit: nSupported
  path "svc_permute_replicate_${x}_gate.csv", emit: gate
  path "svc_permute_replicate_${x}_gate_null.csv", emit: gateNull
  path "svc_permute_replicate_${x}_gate.png", emit: gatePlot

  script:
  """
  # Fit adaptive SVC RFE-CV to one shared phenotype-label permutation
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
    --gate false \
    --gate-permutations '${params.svc_gate_permutations}' \
    --gate-alpha '${params.svc_gate_alpha}' \
    --gate-seed '${params.svc_gate_seed}' \
    --gate-early-stop '${params.svc_gate_early_stop}' \
    --analysis-type permutation \
    --replicate '${x}' \
    --output-prefix 'svc_permute_replicate_${x}'
  """
}
