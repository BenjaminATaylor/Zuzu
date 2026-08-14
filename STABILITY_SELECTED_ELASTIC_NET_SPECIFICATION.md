# Stability-selected elastic-net specification

Status: provisional exploratory integration. Standalone real-data and synthetic pilots and reduced end-to-end smoke tests have been completed, but the method has not yet been validated as a general-purpose feature-selection procedure. The pipeline branch is deliberately opt-in and its settings remain open to later review.

## 1. Purpose

The scientific question for this method is:

> Under which data conditions can a stability-selected elastic-net logistic regression recover phenotype-associated genes reliably, while avoiding unstable selections under weak or absent signal?

This is deliberately different from the question answered by DESeq2, edgeR, and Wilcoxon. Those methods test genes individually or marginally. Elastic-net logistic regression selects genes jointly because they contribute to a sparse additive predictor of phenotype. Its output is therefore a **stable predictive feature signature**, not a conventional DEG list.

Zuzu may nevertheless compare its selected gene list with the lists from the conventional methods, provided that the difference in estimands remains explicit. In the existing full-synthetic benchmark, all methods can be assessed against the same simulator-defined truth set using the same gene universe and confusion-matrix metrics.

This specification does not reactivate, modify, or depend on the SVC branch.

## 2. Current integration scope

The initial pipeline integration has three opt-in components.

| Component | Scope | Existing workflow affected |
| --- | --- | --- |
| Baseline | Run the selector once on the cleaned observed dataset | New elastic-net branch only |
| Constrained null | Rerun the complete selector for every existing shared phenotype-label permutation | New elastic-net branch and aggregate null plots only |
| Full synthetic | Run the selector on every dataset in the existing full-synthetic grid and score its stable signature against the existing truth | New elastic-net branch and aggregate synthetic plots only |

Quasi-permutation remains deferred until its interpretation for a covariance-dependent multivariable selector has been reviewed.

This integration must not:

- alter the numerical inputs, models, thresholds, or results of DESeq2, edgeR, or Wilcoxon;
- alter `DATA_FULLSYNTH` or expand its simulation regimes;
- add a second permutation gate analogous to the SVC gate;
- use the dormant `mlstep` option to enable elastic net;
- add elastic-net work unless the dedicated `enetstep` option is enabled;
- make independent-prediction claims without an independent test set or outer validation procedure.

The existing full-synthetic datasets contain marginal mean-shift truth and are adequate for the initial gene-recovery comparison. They do not assess interaction-only, nonlinear, dispersion-only, or distribution-shape signals.

## 3. Terminology and permitted claims

The following terms must remain distinct:

- **Conventional DEG:** a gene called by a gene-wise test at a stated adjusted-p-value threshold.
- **Base-selected feature:** a gene with a non-zero elastic-net coefficient in one half-sample fit under the fixed base-selector rule.
- **Selection probability:** the observed proportion of half-sample fits in which a gene is base-selected.
- **Stable feature:** a gene whose selection probability meets the prespecified stability cutoff.
- **Stable signature:** the complete set of stable features from one dataset. This set may be empty.
- **Truth gene:** a gene marked as differential by the synthetic generator.

The initial implementation may report:

> Genes repeatedly selected by a sparse additive logistic model under complementary half-sample perturbations.

It must not claim that:

- a selection probability is a p-value, posterior probability, or probability that a gene is biologically relevant;
- the fixed stability cutoff provides BH false-discovery-rate, PFER, or gene-wise error control;
- a selected coefficient is a log fold change or a causal effect;
- an omitted gene lacks marginal differential expression;
- every selected gene is individually necessary for prediction;
- successful phenotype prediction and exact recovery of synthetic truth are equivalent outcomes.

If a signed coefficient is added in a later stage, it must be labelled as a conditional model diagnostic. Correlated genes can change one another's coefficient magnitudes and signs.

## 4. Input and preprocessing contract

### 4.1 Analysis universe

For observed, constrained-null, and any future quasi-permutation analyses, the elastic-net branch will use the same gene universe emitted by `CLEANINPUTS` for the conventional methods. It must not add a phenotype-aware gene filter.

For the existing full-synthetic branch, it will use the complete filtered matrix emitted by `DATA_FULLSYNTH` and the matching saved truth set. The current difference between the observed-data and synthetic filtering paths is outside this specification.

The baseline result table must retain the complete applicable gene universe, including genes that are never selected or are constant within one or more half-samples. Outer permutation and synthetic tasks may retain stable-only tables plus their complete input matrices and compact diagnostics, avoiding thousands of redundant all-gene tables.

### 4.2 Expression transformation

The current transformation is:

1. calculate each sample's library size from the available pre-analysis count universe;
2. calculate counts per million;
3. transform the values as `log2(CPM + 1)`;
4. restrict the matrix to the applicable filtered gene universe;
5. allow `glmnet` to centre and scale predictors using only the samples in the current half-sample fit.

The count columns and phenotype rows must be put in identical sample order and checked explicitly. Reference samples are encoded as 0 and alternative samples as 1. Only an independent two-group design without covariates is initially supported.

Library-size calculation and the CPM transformation are sample-specific and may be performed once before stability resampling. Predictor centring, scaling, zero-variance detection, class weights, and model fitting must be recomputed inside every half-sample.

## 5. Current provisional selection procedure

### 5.1 Provisional exploratory settings

The current permissive procedure uses:

- binomial elastic-net logistic regression fitted with `glmnet`;
- elastic-net mixing parameter `alpha = 0.25`;
- 50 complementary pairs, producing 100 half-sample fits;
- no per-fit feature budget (`q` is unset);
- a fixed inclusive stability cutoff `pi_threshold = 0.90`;
- phenotype-balanced observation weights calculated afresh within every half-sample;
- stratified cross-validation within every half-sample, selecting `lambda.1se` by binomial deviance;
- five requested cross-validation folds, reduced to the largest feasible value for a small half-sample;
- 100 requested lambda values on the `glmnet` path;
- `lambda_min_ratio = 0.001`;
- a fixed and recorded master seed;
- one computational thread per `glmnet` fit, with parallelism applied only at a declared outer level.

These settings deliberately favour a less sparse base selector than the original `alpha = 0.5`, `lambda_min_ratio = 0.01` pilot. They are starting settings for exploratory plots, not validated constants or implicit defaults for future confirmatory work. Section 8 defines how changes must be recorded and evaluated.

`alpha` must satisfy `0 < alpha <= 1`. Ridge regression at `alpha = 0` does not create a sparse selected set and is not part of this procedure.

### 5.2 Complementary half-samples

For each of the 50 pairs, samples will be split without replacement within each phenotype group. For each class, one half receives `floor(n_class / 2)` members and the other receives `ceiling(n_class / 2)`. When one or both class sizes are odd, the seeded allocation of the extra class member is balanced between the two halves so that their total sizes differ by no more than one. Thus the two halves are disjoint, together contain every sample, and remain as close as possible in both total size and phenotype composition.

This phenotype-stratified construction ensures that every fit contains both phenotypes. The implementation must reject a dataset with fewer than six samples in either group, because each half needs at least three members of each class. The grouped inner-cross-validation rule below also requires at least nine samples in every half (at least 18 samples overall), so satisfying the per-group minimum alone is not always sufficient. Paired, repeated-measure, family-related, clustered, or otherwise dependent samples require group-aware resampling and are outside the initial implementation.

All pair allocations must be generated before fitting. The baseline allocations are retained explicitly; compact outer runs retain their labels and seed so the deterministic allocations can be reproduced without publishing a large allocation table for every dataset. Given the same inputs and seed, serial and parallel runs must produce identical allocations and selections.

The randomisation stream must use a canonical phenotype-label order that does not depend on which label is designated as the reference. Reversing the requested contrast must therefore preserve the half-samples, inner folds, selection probabilities, and stable signature; only the outcome encoding and any future coefficient signs may reverse.

The selection count denominator is always 100 for a full exploratory run. A reduced smoke setting uses exactly twice its configured pair count. A failed fit must not be silently removed from the denominator or treated as an empty selection; it makes the analysis fail with a diagnostic identifying the pair and half.

### 5.3 Base elastic-net selector

Each half-sample fit will use:

- `family = "binomial"`;
- `alpha = 0.25` initially;
- an intercept that is never counted as a feature;
- `standardize = TRUE` within that half-sample;
- observation weights of `n_half / (2 * n_class_in_half)`, so the two phenotypes contribute equal total weight;
- a deterministic regularisation path with 100 requested values and `lambda_min_ratio = 0.001`;
- stratified inner cross-validation using binomial deviance, `alignment = "lambda"`, and grouped fold summaries;
- the non-zero coefficients at `lambda.1se` as the base-selected set.

The requested inner fold count is five. For each half-sample, the resolved count is the minimum of the requested count, the number of samples in its smallest phenotype class, and `floor(n_half / 3)`. The final term prevents `cv.glmnet` from automatically replacing grouped fold summaries with ungrouped observation-level summaries when fewer than three observations are available per fold. The existing 11-reference plus 9-alternative synthetic cell consequently uses three folds after complementary splitting. The resolved count must be at least three, every validation fold must contain both phenotypes, and the requested and resolved counts must be retained in diagnostics.

There is no `q` cap. A base fit may therefore select any number of genes supported by its fitted path and `lambda.1se`. This removes a prima facie upper bound on the reported signature, but it does not guarantee a large result: penalisation, correlated substitutes, and instability between half-samples can still make the stable set small.

The implementation will use a small, tested wrapper around `glmnet` so that the binomial family, elastic-net mixing, weights, fold allocation, path settings, selected lambda, and coefficient extraction are explicit. `lambda.min` and path-endpoint behaviour may be retained as diagnostics but do not define the current selected set.

Genes with zero variance in a particular half-sample are unselected for that fit but remain in the full output universe and keep the full denominator of 100 fits. Numerical or convergence failures are errors, not valid empty base sets.

### 5.4 Stable set and current cutoff rule

For gene `j`, define:

```text
selection_probability[j] = selection_count[j] / 100
```

The primary stable set is:

```text
stable_features = genes with selection_probability >= pi_threshold
```

Equality with the threshold is included. The result is valid when this set is empty and must then receive a status such as `ok_empty`, distinct from analytical failure.

The current cutoff is fixed at `pi_threshold = 0.90`; it is not estimated separately for each observed, permuted, or synthetic dataset. Earlier pilots explored null-calibrated cutoffs, but repeating a complete calibration inside every pipeline dataset would multiply the selector workload by approximately 100 and is not part of this provisional integration.

Because `q` is unset, the current procedure does not calculate or claim a theoretical per-family error-rate bound. The stability cutoff is not BH FDR control, a gene-wise significance threshold, or a PFER guarantee. The existing outer constrained-label permutations are instead used as an empirical diagnostic of how often the complete fixed procedure produces stable features under deliberately weakened phenotype alignment.

The metadata must state that the cutoff is fixed, `q` is unset, and no PFER or FDR guarantee is claimed. It must also retain the realised mean, minimum, and maximum base-set sizes so that sparsity and path-boundary behaviour remain visible.

### 5.5 No separate signal gate

Unlike recursive SVC selection, the stability threshold can naturally return an empty signature. The initial method will therefore not add a model-level label-permutation gate around each stability-selection run.

The existing constrained-null branch remains a repeated diagnostic of what the complete fixed selector does after phenotype relabelling. It must not be presented as a gene-wise significance test, a second calibration layer, or an estimate of the fraction of observed calls that are false.

## 6. Output contract

The baseline selector run must retain a per-gene table over the complete analysis universe containing at least:

- gene identifier;
- selection count;
- total half-sample fits;
- selection probability;
- stable-feature status;
- a deterministic rank by decreasing selection probability, with gene identifier used to break ties.

The baseline run must also retain:

- a pair-and-fold allocation table containing pair, half, sample identifier, phenotype, and inner-fold assignment;
- a per-fit diagnostic table containing pair, half, class counts, requested and resolved fold counts, selected count, `lambda.1se`, `lambda.min`, path endpoints, zero-variance count, and completion status;
- run metadata containing the transformation, library-size provenance, class encoding, `alpha`, fixed cutoff, explicit unset `q`, absence of a PFER/FDR claim, pair and fit counts, cross-validation and path settings, weighting rule, random seed, thread settings, and R/package versions;
- a plot of selection probability with the actual cutoff and stable-feature count shown;
- a concise selected-count output for pipeline aggregation.

Every compact outer permutation or synthetic run must instead retain its stable-only gene table, concise summary, fit diagnostics, and metadata. Its input count matrix remains the authoritative tested universe. Pair allocations are deterministic from the retained labels and seed but are written explicitly only for the detailed baseline run. This compact contract prevents output staging from dominating a grid containing hundreds of independently fitted datasets.

The baseline complete per-gene table is its authoritative result. Plots may show only the most highly ranked genes when necessary for readability.

Stability selection produces a selected set, not one uniquely defined final prediction model. No arbitrary half-sample model may be published as “the fitted elastic-net classifier.” A later prediction stage would need a separately specified full-training-data refit and an independent validation design.

## 7. Pipeline integration

### 7.1 Baseline analysis

The baseline branch will run the complete procedure once on the cleaned observed data. Its stable signature is included in the baseline UpSet comparison as `Elastic net stable features`, rather than being labelled as `DEGs`. The exploratory negative-binomial-fit violin diagnostic uses the same four-list membership counts.

The initial integration should publish its own tables, diagnostics, metadata, and plot in a dedicated output directory. It is enabled only by `enetstep` and must remain separately selectable from SVC's `mlstep`.

### 7.2 Constrained-null analysis

For every existing phenotype-label permutation, the complete fixed stability-selection procedure must be rerun using that replicate's shared relabelled dataset. The implementation must not reuse observed-data selection probabilities, class-stratified pair allocations, cross-validation folds, chosen lambdas, or a selected feature set. It does reuse the prespecified `alpha`, path settings, and `pi_threshold = 0.90`, because the outer branch diagnoses that fixed procedure rather than retuning it.

The branch will report at least:

- the number of stable features per replicate;
- the proportion of replicates with an empty stable set;
- each gene's frequency of becoming stable across outer null replicates;
- the distribution of maximum selection probability;
- all resolved settings and seeds.

These outer replicates, not the 100 internal half-samples, are the independent units for summaries of null behaviour.

### 7.3 Existing full-synthetic analysis

The complete selector will be run independently on every existing full-synthetic dataset. The stable-feature vector will be compared with the saved synthetic truth over the full emitted gene universe using the same true-positive, false-positive, true-negative, false-negative, power, FDP, and Youden's-index definitions used elsewhere in Zuzu.

The combined plots and tables must identify the method as a predictive feature selector. Exact truth recovery and any later prediction result are separate outcomes: a correlated proxy can be useful for prediction while still being a false positive against an exact simulator-defined truth set.

No elastic-net setting may be selected by looking at the truth labels of a synthetic replicate later used to report confirmatory performance.

### 7.4 Quasi-permutation is deferred

Quasi-permutation is not part of the first integration. Its current row-wise shuffling destroys gene-gene covariance, while elastic net selects genes conditionally using that covariance. A quasi result could therefore combine feature-selection instability with an intentional change to the method's multivariable input structure.

If the initial integration results justify proceeding, a later quasi-permutation specification must define:

- a method-specific high-confidence anchor rule, such as a prespecified stability threshold stricter than the primary threshold;
- how pre-filter library sizes are handled after quasi-permutation;
- whether destroying background covariance answers a useful self-consistency question for this method;
- how correlated substitute genes are interpreted;
- that the result is method-specific self-consistency, not comparative accuracy.

The anchor threshold must not be chosen retrospectively to obtain a convenient number of anchors.

## 8. Parameter changes and tuning

The current values are explicit so that the exploratory integration evaluates one reproducible procedure. They may be changed when evidence gives a scientific or practical reason to do so.

At minimum, the implementation should expose and validate:

- elastic-net `alpha`;
- complementary-pair count;
- stability cutoff;
- requested inner cross-validation folds;
- class-weighting rule;
- random seed;
- lambda selection rule, path resolution, and lower endpoint;
- allocated CPUs.

The following rules apply:

- `q` remains explicitly unset unless a later specification introduces a separately justified capped procedure;
- changing the cutoff does not create a PFER or FDR guarantee;
- a fixed-parameter run may use either the provisional settings or explicitly supplied settings;
- a tuned run must declare its candidate settings, objective, tie rule, and seeds before fitting;
- the complete tuning rule must be repeated independently in every outer null or synthetic replicate;
- observed agreement with DESeq2/edgeR/Wilcoxon, synthetic truth used for final evaluation, and held-out prediction results must not be used to select the reported settings;
- development simulations may guide settings only when a separate, independently generated confirmatory batch is retained;
- materially different settings are different selection procedures and their metadata and results must not be pooled silently.

The present inner cross-validation selects lambda only. It does not tune `alpha`, `pi_threshold`, the lambda-path bounds, or the number of complementary pairs. Any later tuning of those quantities must be nested inside each outer dataset used for evaluation, or developed on a separate batch and frozen before confirmatory evaluation. Its runtime cost must be reassessed.

## 9. Regression and smoke coverage

The standalone and opt-in smoke tests must verify that:

- count columns and phenotype rows are aligned exactly;
- reference and alternative labels are encoded explicitly;
- the full configured algorithm performs exactly twice the configured number of complementary-pair fits;
- the two halves of every pair are disjoint, complementary, and contain both phenotypes;
- predictor scaling and class weights are fitted within each half-sample;
- inner folds are stratified, are regenerated within each half-sample, and retain grouped fold summaries without exceeding the requested five;
- the 11-reference plus 9-alternative full-synthetic cell resolves to a valid three-fold design where required;
- the selected coefficients come from that fit's `lambda.1se` and `q` remains unset;
- selection probabilities equal selection counts divided by the complete fit count;
- a gene exactly at the stability cutoff is selected;
- zero-variance genes remain in the output universe and are handled consistently;
- empty stable signatures are valid and analytically distinct from failed fits;
- serial and parallel runs give identical results from the same seed;
- swapping reference and alternative labels leaves selection probabilities unchanged within numerical tolerance;
- both directions of phenotype imbalance work;
- unsupported small-class designs fail clearly;
- every outer null replicate reruns the complete selector while retaining the fixed `pi_threshold = 0.90`;
- synthetic confusion counts and metrics agree with a direct calculation;
- no elastic-net work occurs when its dedicated option is disabled;
- inspectable tables, metadata, allocations, diagnostics, and plots are retained.

The ordinary smoke test should continue to omit this branch. An opt-in elastic-net smoke mode may use fewer complementary pairs for structural testing, but its outputs must state clearly that the reduced setting is not the full exploratory procedure.

## 10. Runtime and implementation gate

With 50 complementary pairs, one dataset requires 100 cross-validated `glmnet` fits. At `nperms = 100`, baseline, constrained-null, and the current three-by-three full-synthetic grid together contain 1,001 datasets and therefore approximately 100,100 half-sample fits before any quasi analysis. The test5 pilots suggest roughly 20--30 core-seconds per complete selector, giving a preliminary estimate of approximately 6--8 core-hours for that fixed-cutoff workload. Actual wall time depends on feature count, sample size, executor concurrency, and available CPUs.

Repeating 99 cutoff-calibration relabellings inside every dataset would raise the same workload to approximately 10,010,000 half-sample fits, or roughly 600--800 core-hours by that preliminary estimate. This is why the provisional integration fixes `pi_threshold = 0.90` and uses the existing outer constrained-null branch as a diagnostic rather than nesting another calibration layer.

Runtime assessment must cover:

1. one observed `test5` analysis;
2. a small set of constrained-null analyses;
3. representative full-synthetic datasets at the smallest and largest sample/depth settings;
4. serial and permitted parallel execution.

The standalone pilots have timed observed test5 and representative synthetic analyses. Integrated runs must continue to record timings, particularly for the smallest and largest full-synthetic cells, so the complete `nperms = 100` estimate can be revised. Nested parallelism must be prevented: either Nextflow parallelises independent datasets or a task parallelises its half-sample fits using explicitly allocated CPUs.

The pipeline image must explicitly install `r-cran-glmnet` for fitting and `r-cran-jsonlite` for metadata. The current uncapped fixed-cutoff procedure does not require `stabs`. Exact R and package versions must be recorded in each run's metadata; the repository-wide dependency-pinning limitation remains documented separately.

## 11. Exploratory integration acceptance criteria

The completed standalone pilots must continue to provide evidence that the selector:

- obeys its sampling, cross-validation, fixed-cutoff, empty-result, and reproducibility contracts;
- returns understandable results on observed and clearly null fixtures;
- retains complete diagnostics;
- has a plausible projected runtime.

The opt-in integrated branch is acceptable for exploratory use only if it:

- commonly returns no stable features on clear null data or otherwise has empirically acceptable, fully reported null behaviour;
- recovers at least some planted signals in relevant existing synthetic conditions;
- leaves all non-elastic-net results unchanged;
- remains computationally practical at the intended replicate count;
- produces outputs whose different interpretation from conventional DEGs is clear.

Failure of these criteria should lead to review or abandonment of the method rather than adding an arbitrary top-gene rule or presenting the exploratory cutoff as formal error control.

## 12. References

- Zou H, Hastie T. [Regularization and variable selection via the elastic net](https://academic.oup.com/jrsssb/article/67/2/301/7109482). *Journal of the Royal Statistical Society: Series B*. 2005.
- Friedman J, Hastie T, Tibshirani R. [Regularization Paths for Generalized Linear Models via Coordinate Descent](https://www.jstatsoft.org/article/view/v033i01). *Journal of Statistical Software*. 2010.
- Meinshausen N, Buehlmann P. [Stability selection](https://stat.ethz.ch/Manuscripts/buhlmann/stabilityselection.pdf). *Journal of the Royal Statistical Society: Series B*. 2010.
- Shah RD, Samworth RJ. [Variable selection with error control: another look at stability selection](https://www.statslab.cam.ac.uk/~rds37/papers/Shah%20Samworth%202013%20Variable%20selection%20with%20error%20control%20-%20another%20look%20at%20stability%20selection.pdf). *Journal of the Royal Statistical Society: Series B*. 2013.
- [`glmnet` documentation](https://glmnet.stanford.edu/articles/glmnet.html).
