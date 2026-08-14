process RESAMPLE_SUMMARISE {

  label 'resample_summary'

  publishDir "$params.outdir/subsampling", pattern: '*.csv', mode: 'copy'
  publishDir "$params.outdir/subsampling", pattern: '*.pdf'

  input:
  path summaries, stageAs: 'summaries/summary_?.csv'
  path genes, stageAs: 'genes/genes_?.csv'
  path allocations, stageAs: 'allocations/allocations_?.csv'
  path libraries, stageAs: 'libraries/libraries_?.csv'
  path baseline_membership

  output:
  path "*_subsample_results.csv", emit: results
  path "*_subsample_gene_frequencies.csv", emit: frequencies
  path "sample_subsample_allocations.csv", optional: true, emit: sample_allocations
  path "depth_subsample_library_metrics.csv", optional: true, emit: depth_libraries
  path "*_subsample_plot.pdf", emit: plots

  script:
  """
  # Combine replicate-level outputs and compare them with the full-data lists
  zuzu_resample_summarise.R \
    --summaries summaries \
    --genes genes \
    --allocations allocations \
    --libraries libraries \
    --baseline '${baseline_membership}'
  """
}
