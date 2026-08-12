# Known issues

Zuzu is an in-development research pipeline. Active issues below apply to the current `dev` implementation and are grouped from highest to lowest severity. Resolved issues are retained at the end for traceability.

- **Critical:** can corrupt source data or invalidate a central benchmark.
- **High:** can produce materially incorrect or misleading scientific results.
- **Medium:** limits validity, reproducibility, portability, or output quality under a narrower set of conditions.
- **Low:** incomplete functionality, misleading presentation, or development/documentation debt.

Benchmark results should still be treated as exploratory; the remaining limitations below define where stronger claims are not yet supported.

## Critical

None currently.

## High

None currently.

## Medium

### KI-002: Quasi-permutation self-consistency is labelled as comparative accuracy

The quasi-permutation branch deliberately derives a separate set of high-confidence anchor genes with each method, retains half of those genes unchanged, and asks the same method to recover them after the remaining genes are permuted. This is a method-specific self-consistency analysis, not a benchmark against independent biological truth. Different methods are not expected to share an anchor set or quasi-permuted dataset, and their resulting values are not direct measures of comparative accuracy.

The README now states this distinction, but process outputs still use the fields `power`, `FDP`, and `youden_index`, plots still present Power, FDP, and Youden's index side by side across methods, and the workflow uses `trueDEGs` for method-selected anchors. KI-007 corrected the former ROC AUC terminology, but this presentation can still imply externally validated accuracy in the quasi branch.

Impact: quasi-permutation results can still be mistaken for estimates of biological sensitivity, false discovery, or cross-method performance. Outputs and plots should use anchor-recovery and background-call terminology and make their within-method interpretation explicit.

### KI-013: SVC has no no-signal gate or gene-wise error control

Stage 1 now labels the output as SVC-RFECV candidate features, keeps it distinct from conventional DEG calls, and describes the selected count as the smallest cross-validated maximum rather than an inflection point. However, RFECV still retains at least one feature, its internal maximum is not independent predictive performance, and it produces neither gene-wise p-values nor false-discovery-rate control.

Impact: Stage 1 can report an exploratory predictive candidate set, but it cannot determine that phenotype signal exists or promote that set to a supported signature. The separately specified Stage 2 no-signal gate is still required before that stronger interpretation.

### KI-014: Only simple unpaired two-group designs are supported

DESeq2 and edgeR fit only `~phenotype`, while Wilcoxon performs an unpaired two-sample rank test; all branches require exactly two groups. There is no representation of pairing, donor, batch, repeated measures, interactions, or other covariates.

Impact: the pipeline is conditionally appropriate only for independent two-group experiments; applying it to more complex designs can confound phenotype effects.

### KI-015: The baseline methods do not use an identical testing procedure or estimand

DESeq2 applies its default independent and Cook's filtering. Its default `alpha = 0.1` is used to optimize the independent-filter cutoff even though calls are later made at 0.05; this does not change the stated 0.05 call threshold, but matching `alpha` can improve consistency and power at that threshold. edgeR and Wilcoxon apply BH over all genes retained by the common prefilter. Wilcoxon is an unpaired rank-distribution test on TMM-normalized CPM rather than a negative-binomial mean model; its output contains only adjusted p-values and omits the raw p-value, statistic, effect size, and direction. Its default exact/asymptotic selection can produce extensive tie warnings on tied RNA-seq values, as observed in the current fixture; this is a consistency and reporting concern rather than inherently an invalid test.

Impact: the methods do not have precisely the same multiplicity universe or biological estimand. This can be defensible as a comparison of package-default workflows, but it must be explicit and DESeq2's `alpha` should match the reporting threshold.

### KI-016: The synthetic model only partially reflects the source data

[`DATA_FULLSYNTH`](modules/data_fullsynth.nf) receives the original raw count frame rather than `CLEANINPUTS` output, so its starting gene universe and filtering path differ from the baseline and quasi branches. It resembles the input primarily through its vector of gene row means. Dispersion, effect-size behaviour, and other properties remain compcodeR defaults, while sequencing depths are hard-coded to `1e5`, `1e6`, and `1e7` rather than estimated from the source data.

In the currently resolved compcodeR scalar-effect implementation, differentially expressed variables are assigned from leading variable indices. Zuzu passes means in source row order without randomizing that order, so simulated truth can depend on arbitrary input gene ordering. This detail is also vulnerable to the dependency drift in KI-020.

Impact: results describe this specific compcodeR simulation and preprocessing regime, not necessarily the user's experiment.

### KI-017: Default simulation counts provide little precision and no uncertainty

The default is seven permutations. This is too sparse for stable tail behaviour or calibration near a 0.05 threshold, and synthetic plots report means without intervals or other Monte Carlo uncertainty.

Impact: apparent differences between methods may be simulation noise. Increasing `nperms` helps, but uncertainty should still be reported.

### KI-018: Input and count validation is incomplete

The workflow does not fully enforce file existence, required column names, non-negative finite integer counts, nonempty filtered data, unique sample identifiers, nonblank `sample` and `phenotype` values, or adequate per-group replication. Fractional values are silently rounded. Alignment currently depends on independently sorting sample names and then using positional correspondence; downstream readers also differ in whether they rewrite non-syntactic names.

Impact: malformed or unsuitable inputs can fail late, be silently altered, or produce misaligned or invalid analyses. Estimated counts may also need a method-aware import route rather than unconditional rounding.

### KI-019: Required parameters and value ranges are not fully validated

Only `samplesheet`, `countsframe`, and `reflevel` receive explicit presence checks. The following lack complete defaults or validation:

- `outdir`
- `truecutoff`, including `0 < truecutoff < 1`
- `synthdeg_pct` when `synthstep` is enabled, including a valid proportion range
- `nperms`, including a positive-integer check
- Boolean types for `mlstep`, `synthstep`, and `pareml`; arbitrary nonempty strings can be truthy in workflow conditionals while comparing unequal to `"true"` inside generated scripts

The Stage 1 SVC parameters now have shared defaults and validate `C`, class weighting, scoring, fold count, seed, RFE step, and job count before fitting, but this does not validate the remaining workflow-wide parameters above.

Impact: invalid values can produce malformed scripts, unexpected task counts or seeds, inconsistent branch behaviour, misleading benchmarks, or late task failures.

### KI-020: Runtime dependencies are not reproducibly pinned

The workflow uses the mutable container tag `benjaminataylor/zuzu:latest`. Nextflow is not version-constrained, Python packages have no pinned versions, apt-installed R packages have no fixed versions, and source-installed R packages can resolve newer releases during image builds. SVC prediction bundles record their dependency versions, but their Python pickle format is environment-coupled: a bundle can fail to load under incompatible NumPy or scikit-learn versions and should therefore be used with matching compatible versions.

Impact: rebuilding or rerunning later can change behaviour and numerical results without a source-code change, and saved SVC models are not guaranteed to load in a different software environment.

### KI-022: Automated regression coverage is limited and is not run in CI

[`tests/smoke.sh`](tests/smoke.sh) provides a lightweight default cold integration check for cleaning, QC, the three non-ML methods, permutation, quasi-permutation, aggregation, and plotting. `SMOKE_SYNTH=1` additionally exercises non-ML full-synthetic generation, analysis, aggregation, and plotting. `SMOKE_SVC=1` exercises the baseline and constrained-permutation SVC paths, including their estimator structure, normalization, class direction, outputs, and saved raw-count prediction bundles. The checks assert structural invariants and exact metric relationships, but they do not pin expected inferential results, automate the combined SVC-plus-synthetic run, or run in a CI workflow. The scripts under `sandbox/` remain exploratory rather than runner-integrated tests.

Impact: major workflow breakage and metric inconsistencies are easier to detect, but dependency-driven numerical changes, full-synthetic SVC regressions, and failures outside the covered invariants can still go unnoticed.

### KI-023: Run provenance is incomplete

The workflow does not enable a standard trace, report, timeline, or manifest that records resolved parameters, input hashes, software versions, and the exact container digest alongside results.

Impact: runs are difficult to audit or reproduce after the work directory and console history are gone.

### KI-024: Conventional basic result tables are not published

Baseline DESeq2, edgeR, and Wilcoxon tables are consumed internally but have no `publishDir` or top-level workflow output. They remain only in Nextflow's work cache even though the README presents baseline DEG calling as a standalone use case. Stage 1 now publishes the SVC candidate table and its supporting artifacts under `outdir/svctables`.

Impact: primary inferential outputs are not durable user-facing deliverables.

### KI-025: Published outputs depend on the work cache

With the current default `publishDir` mode, outputs are symlinks into `work/`.

Impact: moving the repository or cleaning the work directory breaks published results, and retaining results requires retaining the much larger task cache.

### KI-026: Collected permutation paths are handled as rendered strings

[`PERMUTE_PLOTS`](modules/permute_plots.nf) accepts collected file paths as `val` inputs and parses their string representation rather than staging a path collection.

Impact: the implementation relies on the whole work root being visible to tasks and is fragile on scratch/remote executors or paths containing spaces and commas.

### KI-027: Small permutation runs publish a non-PDF with a PDF extension

When `nperms < 5`, [`PERMUTE_HISTS`](modules/permute_hists.nf) saves an RData object as `permute_violins.pdf`.

Impact: the published file is not a valid PDF and will fail in document viewers and downstream tooling.

### KI-028: Dummy ML inputs use a nonexistent sentinel path

[`main.nf`](main.nf) defines `dummypath` with `channel.fromPath("0")`. Because `fromPath` does not check existence by default, this literal emits a path even though the file does not exist, so non-ML aggregate processes do launch; the previous claim that the channel was empty was incorrect. `path` consumers nevertheless stage unused broken sentinel symlinks and depend on never reading them.

Impact: executors that materialize inputs may reject the missing source, and any downstream attempt to read the sentinel will fail. Optional inputs or separate process signatures would be clearer.

### KI-029: The cluster profile does not configure a cluster executor

The `cluster` profile changes SVC parallelism and memory only. It does not select an executor, queue, or scheduler.

Impact: the profile still runs locally unless an external site configuration supplies the missing cluster settings.

## Low

### KI-031: Quasi-permutation SVC support is unfinished

The ML conditional in the quasi branch is empty, no SVC quasi analysis is connected, and [`QUASI_PLOTS`](modules/quasi_plots.nf) plots only DESeq2, edgeR, and Wilcoxon results.

### KI-032: Inter-method SVC support is disabled

[`INTER_OBSERVER_BASIC`](modules/inter_observer_basic.nf) accepts an SVC table but hard-codes its ML handling off, so SVC features are absent from reported overlaps.

### KI-033: Legacy SVM modules are disconnected

The R/e1071 SVM modules remain under `modules/`, but their includes and calls are commented out. They contain unresolved feature-selection TODOs and should not be considered supported.

### KI-034: Pipeline stages cannot be selected independently

`synthstep` can disable the synthetic branch and `mlstep` can disable SVC, but QC, baseline, null-permutation, quasi-permutation, and aggregate plotting otherwise run together. There is no supported baseline-only mode.

### KI-036: Documentation and repository status are incomplete

The README provides the scientific concept but not installation instructions, a complete parameter reference, formal input/output schemas, expected resource use, troubleshooting guidance, or interpretation cautions. Prototype scripts and generated research artifacts are not clearly distinguished from supported pipeline components.

## Resolved

### KI-001: Input staging could overwrite source files

Resolved 2026-08-07. [`CLEANINPUTS`](modules/cleaninputs.nf), [`CREATE_BREAKS`](modules/create_breaks.nf), and [`DATA_FULLSYNTH`](modules/data_fullsynth.nf) now stage direct user inputs under controlled names in a `staged_inputs/` subdirectory. Their fixed task outputs therefore cannot share a path with, or write through, a staged source symlink.

Regression coverage deliberately uses the colliding input names `samplesheet.csv` and `countsframe_clean.csv`, verifies immutable snapshots after the run, and passed under Nextflow's symlink staging mode. An extended synthetic run also passed with inputs named `breaks.csv` and `synthcounts.csv`; the source checksums remained unchanged.

### KI-003: Quasi-permutation sample sizes were not applied downstream

Resolved 2026-08-10. The DESeq2, edgeR, and Wilcoxon quasi-data processes now permute and restore counts only within the selected sample subset and emit the matching subset sample sheet. Their analyzers retain the replicate and both requested group sizes, and fail if those sizes, the observed phenotype counts, the sample-sheet rows, or the count-matrix columns disagree.

The smoke regression uses an intentionally unbalanced 30-reference/24-alternative fixture, exercises all three generated size pairs for all three methods, and verifies every quasi result records a consistent total. The cold end-to-end run passed.

### KI-004: The requested reference level was lost before DESeq2 and edgeR ran

Resolved 2026-08-10. Every DESeq2 and edgeR model boundary now reconstructs the two-level phenotype factor after reading CSV, asserts that `params.reflevel` is the reference, and derives the sole alternative level. DESeq2 result extraction uses an explicit alternative-versus-reference contrast; edgeR tests an explicit numeric contrast after verifying that its phenotype coefficient encodes the same direction. The DESeq2 QC model applies the same reference-level contract.

The smoke fixture now deliberately gives the requested reference a name that sorts after the alternative. Its cold non-ML run passed all basic, permutation, and quasi model assertions, and verifies that a known alternative-enriched gene has positive DESeq2 and edgeR fold changes. A separate reference-swap check confirmed that both methods reverse the effect sign while retaining the same FDR-selected gene set on the fixture.

### KI-005: The null permutation destroyed sample-level experimental structure

Resolved 2026-08-10. [`DATA_PERMUTE`](modules/data_permute.nf) now leaves the cleaned count matrix unchanged and permutes the phenotype labels once per replicate while preserving both group sizes. The resulting keyed label assignment is shared by DESeq2, edgeR, Wilcoxon, and SVC, so every method sees the same negative-control dataset. Each replicate records its seed, labels, group sizes, overlap, and label correlation in `permutation_metadata.csv`.

Assignments are sampled directly and uniformly from the fixed-margin label vectors whose absolute correlation with the original labels is strictly less than 0.6. This preserves the rule against nearly reproducing the original grouping or its reversal without a retry loop. It also handles unbalanced groups, resolving KI-009, and fails immediately with a clear message for the infeasible two-sample 1:1 case. This remains a constrained diagnostic under phenotype-label exchangeability, not unrestricted classical permutation inference, and its discovery counts do not estimate the fraction of baseline calls that are false.

The smoke regression verifies that the forwarded count matrix is byte-identical to the cleaned input, group sizes and sample order are preserved, every label correlation satisfies the constraint, no replacement count matrix is created, and all non-ML methods consume each replicate exactly once.

### KI-006: The common gene filter was phenotype-aware and raw-depth-dependent

Resolved 2026-08-10. [`CLEANINPUTS`](modules/cleaninputs.nf) now calculates counts per million from each sample's rounded, unfiltered library total and retains genes with at least 1 CPM in at least as many samples as the smaller phenotype group. Group size determines the required prevalence, but group membership is not used in any gene's filtering decision. The retained count values are not transformed, and the same filtered matrix is shared by QC and every baseline, null-permutation, and quasi-permutation method. If more than half of the input genes are removed, the process reports the counts and percentage in a visible warning while continuing with the retained genes. The full-synthetic branch remains separate under KI-016.

The smoke fixture includes sentinel genes confirming that below-CPM and insufficiently prevalent genes are removed, sufficiently prevalent genes are retained, expression confined to one complete phenotype group is not discarded, retained counts remain unchanged, and a greater-than-50% removal warning reaches the console. A separate fixed-margin relabelling check confirmed that phenotype assignments do not alter the retained gene set.

### KI-007: ROC AUC was calculated from binary DEG calls

Resolved 2026-08-10. Every quasi and full-synthetic analyzer now calculates signed Youden's index directly from the binary truth and call vectors as sensitivity minus false-positive rate. Outputs use the `youden_index` field, and plots label the metric as Youden's index over its full -1 to 1 range. The pROC calculation and its automatic score-direction selection were removed, so an inverted result remains negative rather than being presented as successful discrimination.

The metric is explicitly documented as a summary at the fixed DEG or feature-selection threshold. Zuzu no longer describes it as a threshold-independent ROC AUC or as a balance between power and false discovery proportion.

### KI-008: Empty discovery or truth sets produced invalid metrics and summaries

Resolved 2026-08-11. Whenever a quasi or full-synthetic analyzer evaluates metrics, it now derives power, FDP, and Youden's index from explicit true-positive, false-positive, true-negative, and false-negative counts over one fixed gene universe. Outputs also retain the tested, called, truth, and background sizes plus a metric status. `n_tested` denotes this full input universe, so genes with unavailable adjusted p-values after internal method filtering remain uncalled rather than leaving the denominator. KI-010 subsequently removed the quasi branch's temporary ten-anchor censoring, so active quasi rows now retain complete call and confusion counts even when few or no anchors are available. A method making no calls receives FDP 0 under the standard `FP / max(calls, 1)` convention; power and Youden's index are also zero when both truth classes exist, but remain undefined when their required classes are absent.

[`DATA_FULLSYNTH`](modules/data_fullsynth.nf) now rejects empty generated matrices, duplicate identifiers, truth genes absent from the generated count matrix, and synthetic benchmarks lacking either a true or background class. The analyzers repeat the truth-universe checks defensively. SVC pare mode records unmodelled genes as uncalled so debugging does not change the benchmark denominator.

[`QUASI_PLOTS`](modules/quasi_plots.nf) publishes the combined per-replicate `quasi_results.csv`. [`FULLSYNTH_PLOTS`](modules/fullsynth_plots.nf) publishes `synth_results.csv` and an untransformed, full-precision `synth_summary.csv`; each summary records `n_valid`, `n_total`, and `n_no_calls`, while plot labels show the valid denominator. Missing metrics are averaged only over explicitly counted finite replicates, with an all-missing group retained as `NA`. The standard smoke test verifies the quasi confusion-matrix identities, metric formulas and bounds, zero- and low-anchor statuses, row keys, and combined output; the opt-in synthetic path verifies the corresponding per-replicate results, shared truth universes, summaries, and plot.

### KI-009: The null-permutation acceptance rule assumed approximately balanced groups

Resolved 2026-08-10 alongside KI-005. The raw 20--80% agreement rule and its 100-attempt retry loop were replaced by exact enumeration of the correlations attainable at the observed binary group sizes. The process draws directly from all valid fixed-margin assignments, so heavily unbalanced designs no longer fail merely because their unavoidable raw agreement is high.

### KI-010: Quasi results were censored when too few retained genes were available

Resolved 2026-08-11. The active DESeq2, edgeR, and Wilcoxon quasi analyzers now run for every retained-anchor count. Rows with one to ten anchors retain their complete confusion counts and calculated metrics with a `low_anchor_count` status when background genes remain, rather than disappearing behind an arbitrary cutoff. When no anchors are retained, the analyzers still report calls, confusion counts, and FDP under a `no_anchors` status; power and Youden's index remain `NA` because there is no positive class. If no background genes remain, Youden's index is `NA` under `youden_undefined_no_background` regardless of anchor count.

[`QUASI_PLOTS`](modules/quasi_plots.nf) distinguishes adequate, low, and empty anchor sets by point shape and states why power and Youden's index are absent for empty sets. The smoke validator enforces the count identities, formulas, bounds, and status precedence for zero-, low-, and adequate-anchor rows, rejects defined recovery metrics without anchors, and guards the active analyzers against reintroducing the old gate.

### KI-011: Synthetic analysis ignored the alternative-group sample size

Resolved 2026-08-11. [`DATA_FULLSYNTH`](modules/data_fullsynth.nf) now asks compcodeR to generate both conditions at the larger requested group size, then uses the replicate seed to select the exact requested number of reference and alternative samples. It verifies both group counts and sample-sheet/count-column alignment before emitting the dataset. The analyzers retain and validate `refnum`, `altnum`, and `totalnum`, while [`FULLSYNTH_PLOTS`](modules/fullsynth_plots.nf) groups results and labels sample sizes using both group counts.

The opt-in synthetic smoke validator checks every generated sample-sheet/count-matrix pair directly. It requires the exact 11+9, 20+16, and 30+24 fixture designs at every depth and replicate, unique sample identifiers, aligned count headers, the expected task and pair frequencies, and consistent sample-size fields in the per-replicate and summary tables. The cold non-ML synthetic smoke run passed.

### KI-012: SVC preprocessing leaked information across validation folds

Resolved 2026-08-11. The three active SVC modules now call one shared implementation that applies `log2(CPM + 1)` and passes unscaled expression into `RFECV`. Its estimator is a scikit-learn `Pipeline` containing `StandardScaler` followed by a linear SVC, so every feature-elimination fit learns scaling from its own training fold. Ordinary and constrained-permutation fits use the pre-filter library sizes already calculated by `CLEANINPUTS`; full-synthetic fits use all genes in the unchanged emitted synthetic count matrix.

The procedure explicitly maps the requested reference to class 0, uses balanced class weights and balanced accuracy, and records shuffled seeded stratified folds. It rejects a phenotype group smaller than the requested fold count rather than silently constructing an unsupported split. Its fixed defaults are exposed as validated parameters and retained in metadata rather than treated as immutable constants.

Every run now publishes a complete candidate-feature table, fold assignments, curve data, plot, metadata, selected count, and a model bundle containing the fitted scaler-plus-SVC and the gene universe needed to reproduce raw-count normalization. The opt-in SVC smoke passed under the pipeline container and verified both imbalance directions, the small-class failure, and raw-count prediction from the persisted bundle.

### KI-021: Current HEAD lacked a documented clean end-to-end run

Resolved 2026-08-11. A cold Stage 1 run completed the baseline, constrained-permutation, quasi-permutation, and existing full-synthetic branches with `mlstep=true`, `synthstep=true`, `pareml=true`, and one replicate. It executed one baseline SVC fit, one constrained-permutation fit, all nine SVC full-synthetic fits, and all nine SVC postprocessors without cache reuse. The combined synthetic table contained nine rows for each of DESeq2, edgeR, Wilcoxon, and `SVC-RFECV candidates` over the unchanged 11+9, 20+16, and 30+24 by three-depth grid.

This establishes mechanical completion of the current pared workflow, not benchmark-scale runtime, numerical reproducibility across dependency versions, or scientific validation of SVC as a DEG procedure.

### KI-030: SVC feature-selection plots used the wrong x-axis values

Resolved 2026-08-11. SVC curve tables and plots now use the actual feature counts evaluated by RFECV. The implementation reads scikit-learn's `n_features` result when available and reconstructs the same ascending count sequence for older releases, then verifies that it aligns with the score vector. For example, 100 genes at a step of 50 are plotted at 1, 50, and 100 rather than at evaluation indices 1, 2, and 3.

The selected point is marked explicitly and verified to be the smallest evaluated feature set at the maximum mean fold score. Every fold score and its mean and standard deviation are retained alongside the plot. The implementation passed direct checks under scikit-learn 1.4 and the pipeline image's scikit-learn 1.8.

### KI-035: The quasi plot had an incorrect y-axis label

Resolved 2026-08-10 alongside KI-003. [`QUASI_PLOTS`](modules/quasi_plots.nf) now labels the y-axis as a self-consistency metric and identifies each facet using both the reference- and alternative-group sample counts.
