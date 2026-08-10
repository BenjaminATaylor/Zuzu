# Known issues

Zuzu is an in-development research pipeline. Active issues below apply to the current `dev` implementation and are grouped from highest to lowest severity. Resolved issues are retained at the end for traceability.

- **Critical:** can corrupt source data or invalidate a central benchmark.
- **High:** can produce materially incorrect or misleading scientific results.
- **Medium:** limits validity, reproducibility, portability, or output quality under a narrower set of conditions.
- **Low:** incomplete functionality, misleading presentation, or development/documentation debt.

Until the critical and high-severity issues are resolved, benchmark results should be treated as exploratory.

## Critical

None currently.

## High

### KI-006: The common gene filter is phenotype-aware and raw-depth-dependent

[`CLEANINPUTS`](modules/cleaninputs.nf) retains a gene when the maximum of its two raw phenotype-group means exceeds five. This specific rule is contrast-dependent and ignores library size.

Impact: sequencing-depth differences and chance group imbalances can affect which hypotheses enter testing. A label-blind CPM/prevalence rule, or a clearly documented purpose-built design-aware rule such as edgeR's `filterByExpr`, would provide a more defensible testing universe.

### KI-007: Reported ROC AUC is calculated from binary DEG calls

The quasi and full-synthetic analyzers pass binary selected/not-selected calls to the ROC calculation instead of a continuous statistic or score. With a binary predictor, the resulting AUC is equivalent to balanced accuracy at one threshold, potentially after pROC automatically reverses the score direction. `normAUC = 2 * AUC - 1` is the corresponding single-threshold Youden index rather than an area over a multi-threshold ROC curve.

Affected modules include the `*_quasi.nf`, `*_fullsynth.nf`, and [`svc_fullsynth_postprocess.nf`](modules/svc_fullsynth_postprocess.nf) analyzers.

Impact: `normAUC` is mislabelled and does not provide the claimed threshold-independent balance of power and false discovery.

### KI-008: Empty discovery or truth sets produce invalid metrics and summaries

The quasi and synthetic analyzers calculate FDP as `1 - TP / calls`. When a method makes no calls, this becomes `0 / 0` and yields `NaN`. [`FULLSYNTH_PLOTS`](modules/fullsynth_plots.nf) then averages values without removing missing values, so one empty-call replicate can blank an entire summary cell.

If the simulated truth set is empty, power also divides by zero and ROC has only one response class. This can occur with an invalid zero DEG proportion or if simulator filtering removes all simulated DE genes, and may abort analysis rather than merely produce a missing summary.

Impact: valid empty-set outcomes can corrupt summaries or fail tasks. Counts of TP, FP, TN, FN, tested genes, called genes, and truth genes should be retained; zero-call and zero-truth conventions should be defined explicitly.

### KI-010: Quasi results are censored when too few retained genes are available

The quasi analyzers emit only `NA` metrics when the retained truth set contains ten or fewer genes.

Impact: methods producing fewer high-confidence anchor calls, stringent cutoffs, and smaller sample sizes disappear from summaries rather than recording low or uncertain performance, creating method- and sample-size-dependent missingness.

### KI-011: Synthetic analysis ignores the alternative-group sample size

[`DATA_FULLSYNTH`](modules/data_fullsynth.nf) accepts both `refnum` and `altnum`, but passes only `refnum` to compcodeR as `samples.per.cond`. Every synthetic dataset is therefore balanced even when the source groups and calculated breaks are unequal.

Impact: synthetic benchmarks do not represent the requested unbalanced design and may contain a different total sample count from the source experiment.

### KI-012: SVC preprocessing leaks information across validation folds

The SVC modules fit `StandardScaler` to the complete dataset before `RFECV` performs cross-validation.

Affected modules: [`svc_basic.nf`](modules/svc_basic.nf), [`svc_permute.nf`](modules/svc_permute.nf), and [`svc_fullsynth.nf`](modules/svc_fullsynth.nf).

Impact: validation folds influence their own preprocessing, which can produce optimistic feature-selection scores. Scaling and classification should be composed in a scikit-learn `Pipeline` fitted independently within each fold.

### KI-013: SVC-selected features are interpreted more strongly than the method supports

SVC/RFECV identifies a subset that optimizes cross-validated classification accuracy. It does not identify a statistical inflection point, prove that removed genes contribute only to overfitting, produce differential-expression p-values, or control false discovery rate. RFECV also always selects at least one feature.

Impact: selected genes are predictive features, not directly equivalent to FDR-controlled DESeq2, edgeR, or Wilcoxon DEGs. Comparisons and labels that treat them as the same quantity are misleading.

## Medium

### KI-002: Quasi-permutation self-consistency is labelled as comparative accuracy

The quasi-permutation branch deliberately derives a separate set of high-confidence anchor genes with each method, retains half of those genes unchanged, and asks the same method to recover them after the remaining genes are permuted. This is a method-specific self-consistency analysis, not a benchmark against independent biological truth. Different methods are not expected to share an anchor set or quasi-permuted dataset, and their resulting values are not direct measures of comparative accuracy.

The README now states this distinction, but process outputs and plots still label the quantities as `power`, `FDP`, and `normAUC`, use `trueDEGs` for method-selected anchors, and display methods side by side. The separate mathematical problem with `normAUC` is covered by KI-007.

Impact: quasi-permutation results can still be mistaken for estimates of biological sensitivity, false discovery, or cross-method performance. Outputs and plots should use anchor-recovery and background-call terminology and make their within-method interpretation explicit.

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
- `stepsize` when SVC is enabled
- `njobs`, which is set only by the `local` and `cluster` profiles
- `nperms`, including a positive-integer check
- Boolean types for `mlstep`, `synthstep`, and `pareml`; arbitrary nonempty strings can be truthy in workflow conditionals while comparing unequal to `"true"` inside generated scripts

Impact: invalid values can produce malformed scripts, unexpected task counts or seeds, inconsistent branch behaviour, misleading benchmarks, or late task failures. The default ML configuration can generate invalid Python for `n_jobs` when no profile supplies `njobs`.

### KI-020: Runtime dependencies are not reproducibly pinned

The workflow uses the mutable container tag `benjaminataylor/zuzu:latest`. Nextflow is not version-constrained, Python packages have no pinned versions, apt-installed R packages have no fixed versions, and source-installed R packages can resolve newer releases during image builds.

Impact: rebuilding or rerunning later can change behaviour and numerical results without a source-code change.

### KI-021: Current HEAD lacks a documented clean end-to-end run

The current commit now has a successful cold smoke run covering the non-ML, non-synthetic workflow. The only full default-configuration run used `-resume`, predates the latest plotting/profile changes, and reused more cached tasks than it executed.

Impact: the SVC and full-synthetic branches have not yet been shown to complete together from a cold cache on the current commit.

### KI-022: Automated regression coverage is limited and is not run in CI

[`tests/smoke.sh`](tests/smoke.sh) provides a lightweight cold integration check for cleaning, QC, the three non-ML methods, permutation, quasi-permutation, aggregation, and plotting. It deliberately omits SVC and full-synthetic benchmarking, does not assert numerical scientific results, and is not run by a CI workflow. The scripts under `sandbox/` remain exploratory rather than runner-integrated tests.

Impact: major non-ML workflow breakage is now easier to detect, but regressions in omitted branches or scientific calculations can still go unnoticed.

### KI-023: Run provenance is incomplete

The workflow does not enable a standard trace, report, timeline, or manifest that records resolved parameters, input hashes, software versions, and the exact container digest alongside results.

Impact: runs are difficult to audit or reproduce after the work directory and console history are gone.

### KI-024: Basic result tables are not published

Baseline DESeq2, edgeR, Wilcoxon, and SVC tables are consumed internally but have no `publishDir` or top-level workflow output. They remain only in Nextflow's work cache even though the README presents baseline DEG calling as a standalone use case.

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

### KI-030: SVC feature-selection plots use the wrong x-axis values

The SVC modules label the x-axis as the number of selected features but plot `1..n_scores`. When `RFECV.step` is greater than one, these are evaluation indices rather than the actual feature counts.

Impact: the plotted curve cannot be read as feature count versus validation accuracy and may misstate the location of the selected optimum.

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

### KI-009: The null-permutation acceptance rule assumed approximately balanced groups

Resolved 2026-08-10 alongside KI-005. The raw 20--80% agreement rule and its 100-attempt retry loop were replaced by exact enumeration of the correlations attainable at the observed binary group sizes. The process draws directly from all valid fixed-margin assignments, so heavily unbalanced designs no longer fail merely because their unavoidable raw agreement is high.

### KI-035: The quasi plot had an incorrect y-axis label

Resolved 2026-08-10 alongside KI-003. [`QUASI_PLOTS`](modules/quasi_plots.nf) now labels the y-axis as a self-consistency metric and identifies each facet using both the reference- and alternative-group sample counts.
