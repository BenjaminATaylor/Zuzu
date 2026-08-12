# Staged SVC development and benchmark specification

Status: Stage 1 implemented. Stage 2 remains a separately approved exploratory extension, and Stage 3 remains a separate confirmatory-benchmark proposal.

## 1. Purpose

The eventual scientific question is:

> Under which data-generating conditions does a fully specified SVC-based gene-selection procedure recover simulator-defined differentially expressed genes as reliably as, or more reliably than, DESeq2, edgeR, and Wilcoxon?

Answering that question requires more than correcting the current SVC code. It requires a separate confirmatory benchmark. These two pieces of work must not be implemented as one change.

This document therefore defines three stages:

1. rehabilitate and re-enable SVC within the existing pipeline;
2. add an SVC-specific no-signal gate and perform an initial evaluation;
3. only if justified, implement an expanded confirmatory benchmark.

Each stage ends with a review before the next begins. Stage 1 is the immediate implementation scope. The expanded benchmark in Stage 3 is a future research project, not part of the initial SVC fix.

## 2. Stage boundaries

| Stage | Purpose | Synthetic data | Effect on non-SVC workflow |
| --- | --- | --- | --- |
| 1. SVC rehabilitation | Make the existing classifier mechanically and statistically coherent | Use the existing `DATA_FULLSYNTH` datasets unchanged | None intended |
| 2. Signal gate and initial evaluation | Allow SVC to report no supported signature and assess its behaviour in the current setting | Continue using the existing full-synthetic benchmark, plus small SVC-specific null checks | None intended |
| 3. Expanded benchmark | Test whether SVC can empirically behave as a DEG method across defined conditions | Add independent test cohorts, additional signal regimes, matching nulls, and stability experiments | Conventional methods are rerun inside the new benchmark, but their ordinary workflows remain unchanged |

Through Stages 1 and 2, the following are explicitly out of scope:

- redesigning `DATA_FULLSYNTH`;
- changing the DESeq2, edgeR, or Wilcoxon models or thresholds;
- changing the existing constrained-permutation procedure;
- changing the quasi-permutation procedure;
- increasing the shared `nperms` parameter to benchmark-scale replicate counts;
- replacing existing non-SVC plots or result tables;
- making confirmatory claims that SVC is a valid DEG alternative.

The expensive Stage 3 work must use its own opt-in workflow, parameters, replicate channels, and output directory. It must not expand the runtime of the ordinary pipeline or reuse `nperms` as its simulation count.

## 3. Terminology and permitted claims

The following terms must be kept distinct:

- **Selected gene:** the general term for any gene returned by a method.
- **Conventional DEG:** a gene selected by a gene-wise statistical test at a stated adjusted-p-value threshold.
- **SVC candidate feature:** a gene retained by the RFECV procedure.
- **Supported SVC signature:** the candidate feature set when the Stage 2 model-level signal gate passes.
- **Truth gene:** for the primary synthetic comparison, a gene with a non-zero phenotype effect on its expected normalised abundance.
- **Predictive signature:** a set of genes used together to predict phenotype.

Stage 1 may report:

> SVC-RFECV candidate features selected for discriminating the two phenotype labels.

It may also report exploratory overlap, power, FDP, and Youden's index against the truth supplied by the existing full-synthetic step, provided those results are explicitly described as recovery by a predictive feature selector.

Stage 2 may distinguish an unsupported candidate set from a supported signature. Passing the gate establishes evidence of model-level phenotype structure under its permutation assumptions; it does not give individual genes p-values or FDR control.

Only Stage 3 can support a condition-specific claim that SVC behaved as an empirically DEG-compatible procedure. No stage may support:

- universal error control outside the evaluated regimes;
- gene-level p-values or adjusted p-values for SVC features;
- interpreting an SVC coefficient as a log fold change;
- calling every retained feature individually necessary for prediction;
- calling every removed gene non-predictive or a cause of overfitting;
- treating the internal RFECV score as independent predictive performance;
- claiming that good phenotype prediction proves good DEG recovery.

SVC outputs must use a field such as `selected_feature`, not `DEGstatus`. “More reliable” must always name the relevant outcome: truth recovery, false selection, stability, or independent prediction.

## 4. Stage 1: SVC rehabilitation

### 4.1 Workflow boundary

Stage 1 changes the SVC modules and the code that reads their outputs. It must not alter the numerical inputs or analysis behaviour of DESeq2, edgeR, or Wilcoxon.

The existing `CLEANINPUTS` filtered-count output remains the common analysis universe. SVC additionally needs library sizes calculated before this filter. That information must be supplied through an SVC-specific preparation step or an additional named output without changing the tuple currently consumed by the non-SVC modules.

For the existing full-synthetic branch, the count matrix emitted by `DATA_FULLSYNTH` remains the Stage 1 input universe. Library sizes are calculated from all genes in that emitted matrix. The generator, truth set, sample-size grid, and post-generation total-count filter remain unchanged.

### 4.2 Input transformation

SVC must not be fitted directly to unnormalised raw counts. The rehabilitated procedure will:

1. calculate sample library sizes from all available pre-analysis genes;
2. transform expression to `log2(CPM + 1)`;
3. restrict the matrix to the common filtered gene universe;
4. fit `StandardScaler` separately inside every RFECV training fold.

Scaling and the linear SVC must be contained in a scikit-learn `Pipeline`. Validation samples must not influence their own fitted centring or scaling.

### 4.3 Classifier and recursive elimination

The proposed initial procedure is:

- a linear SVC;
- an RFECV importance getter that reads `coef_` from the fitted SVC pipeline step;
- `C = 1`;
- `class_weight = "balanced"`;
- balanced accuracy as the RFECV score;
- five stratified, shuffled, and seeded folds;
- explicit reference encoding as class 0 and alternative encoding as class 1;
- `step = 50` genes;
- `min_features_to_select = 1`.

Five folds are used only when the smaller phenotype group contains enough independent samples. Unsupported designs must fail clearly. Paired, repeated, or otherwise related samples require group-aware splitting and are outside the initial two-independent-group implementation.

The selected feature count is the evaluated count with the highest mean fold score. Ties are resolved in favour of the smaller set. This is the **CV-optimal feature count**, not a statistical inflection point and not evidence that only overfitting genes have been removed.

Plots must use the actual evaluated feature counts on the x-axis. The elimination step, fold allocation, score, chosen count, and random seeds must be retained.

These are reference defaults for Stage 1, not immutable scientific constants. The implementation should expose `C`, fold count, class weighting, RFE step, and scoring rule as validated parameters and record the resolved values in every run.

Parameters may be changed or tuned when there is a scientific or practical reason. The constraint is on **how** this is done:

- a fixed-parameter run may use either the defaults or explicitly supplied values;
- a tuned run must declare its search space, selection rule, and random seeds before fitting;
- tuning must use training data only and must be nested inside any cross-validation used to estimate prediction;
- the complete tuning and feature-selection search must be repeated for every Stage 2 gate permutation;
- synthetic truth labels, held-out test performance, and confirmatory results must not be used to choose the settings later reported on those same data;
- results produced by materially different settings or tuning rules must retain those configurations and be compared as distinct procedures.

`C` and class weighting are model hyperparameters and can be tuned within that framework. Fold count and scoring metric primarily define the validation question and should normally be chosen from the study design rather than selected because they give the most favourable result. RFE step is principally a runtime-versus-resolution choice, although it can change the selected set and must therefore be recorded.

Stage 1 uses the reference defaults without a hyperparameter search so that the initial change tests the corrected mechanics in isolation. A later tunable mode does not require the expanded Stage 3 benchmark, but it does require the safeguards above.

### 4.4 Stage 1 outputs

Each ordinary SVC run must retain:

- a complete table over the input gene universe with `selected_feature` status;
- selected gene identifiers and final ranks;
- coefficient signs as model diagnostics only;
- actual evaluated feature counts and fold scores;
- the feature-selection plot;
- preprocessing, classifier, CV, and seed parameters;
- the complete fitted transformation-plus-model object needed to score new raw-count samples.

Coefficient signs may change in the presence of correlated predictors and must not be presented as biological effect directions.

### 4.5 Existing synthetic evaluation

Stage 1 passes the rehabilitated SVC through the existing full-synthetic datasets. It does not change `DATA_FULLSYNTH` or rerun the conventional methods under a new simulation design.

The existing power, FDP, and Youden calculations may be applied to the SVC candidate set over the existing emitted gene universe. These are useful exploratory diagnostics, but they do not establish null error control, independent prediction, selection stability, or general validity as a DEG method.

The conventional methods' existing synthetic values should remain unchanged apart from harmless presentation changes required to add or relabel the SVC column.

### 4.6 Stage 1 regression coverage

An opt-in SVC smoke path must verify that:

- count columns and phenotype rows are identically ordered;
- fold-level scaling is not fitted on validation samples;
- actual feature counts are plotted;
- reference and alternative classes are encoded explicitly;
- both imbalanced phenotype directions work;
- unsupported small-class designs fail clearly;
- pare mode preserves the complete output gene universe;
- saved models contain every transformation required for prediction;
- inspectable SVC tables, model metadata, and plots are retained.

The current default smoke test continues to use `--mlstep false` and must remain unchanged in scope and runtime.

## 5. Stage 2: no-signal gate and initial evaluation

### 5.1 Purpose

RFECV always retains at least its minimum feature count. A candidate set must therefore not automatically become a supported result.

Stage 2 adds a model-level gate entirely within the SVC branch. It does not redesign the full-synthetic generator and it does not replace the existing constrained-permutation analysis.

### 5.2 Gate procedure

The proposed gate statistic is the maximum mean fold-level balanced accuracy on the RFECV curve. This optimised score is not an estimate of future predictive accuracy, but it can be used as a permutation-test statistic when the same optimisation is repeated for every permuted label vector.

Uniformly permute the training labels in a way permitted by the study design and repeat the complete SVC procedure, including any label-dependent filtering, scaling, hyperparameter tuning, recursive elimination, and maximisation over feature counts:

```text
p = (1 + number of null scores at least as large as the observed score) /
    (1 + number of null permutations)
```

The proposed initial count is 99 permutations, giving p-value resolution of `0.01`. Observed and permuted runs use the same deterministic split-generation rule and seeds, with stratified folds regenerated for each label vector. A failing gate run may stop once enough exceedances make `p <= 0.05` impossible, but it must then report only `p > 0.05`.

If the gate passes, the complete-data candidate set becomes the supported SVC signature. If it fails, the supported set is empty and the public status is `no_predictive_signal`; the forced candidate set remains available only for diagnostics.

This is a model-level test. It does not attach significance or error guarantees to individual retained genes.

### 5.3 Relationship to the existing permutation workflow

The current constrained permutations remain a self-consistency diagnostic shared with the conventional methods. They must not be used as the gate null because they deliberately exclude label assignments too similar to the observed grouping.

Gate permutations use a separate SVC-specific parameter and channel. The existing SVC constrained-permutation branch may continue to summarize candidate-feature selection, clearly labelled as such; it must not recursively run the full signal gate for every diagnostic permutation.

### 5.4 Stage 2 evaluation and stopping point

Stage 2 uses:

- the ordinary input dataset;
- the existing full-synthetic datasets;
- a small SVC-specific complete-null fixture or targeted null simulation;
- an opt-in smoke configuration with fewer permutations for structural testing only.

The initial review asks:

- Does the corrected SVC complete reliably?
- Does the gate return an empty set on clearly null data?
- Does it retain signal in at least some existing full-synthetic conditions?
- Are candidate and supported outputs understandable and reproducible?
- Is the runtime compatible with a larger study?

The project stops here for review. Stage 3 is undertaken only if the initial behaviour is scientifically interesting and the computational cost is acceptable.

## 6. Stage 3: expanded confirmatory benchmark

Stage 3 is a separate, opt-in research workflow. It does not replace the ordinary baseline, constrained-permutation, quasi, or current full-synthetic branches unless a later explicit decision is made to retire them.

Development and confirmatory simulation batches must be separate. Pilot results may be used to lock procedures and estimate runtime, but must not be included in the final confirmatory estimates.

### 6.1 Truth and evaluation universes

For scenario `theta` and independently generated replicate `r`, define:

- `U(theta, r)`: the complete simulated gene universe;
- `T(theta, r)`: truth genes with a non-zero phenotype effect on expected normalised abundance;
- `S(method, theta, r)`: the genes selected by a locked method.

The primary target is marginal differential expression, matching the main target of DESeq2, edgeR, and Wilcoxon. Dispersion-only, distribution-shape, interaction, or jointly predictive effects without a marginal mean difference require separately named truth sets and cannot be counted as recovery of the primary DEG truth.

Report two universes:

1. **End-to-end universe:** all simulated genes. A truth gene removed by preprocessing is a false negative.
2. **Conditional universe:** genes surviving the common training-derived filter.

The end-to-end result is primary. The conditional result distinguishes filtering losses from selector losses.

### 6.2 Validation design

One independently generated training/test pair is one Monte Carlo replicate. All methods receive the same latent parameters, training samples, phenotype labels, gene universe, truth set, and test cohort.

For each replicate:

1. draw gene means, dispersions, truth identities, and phenotype effects once;
2. generate disjoint training and test samples from those fixed parameters;
3. perform filtering, fitting, selection, tuning, and the SVC gate using training samples only;
4. calculate truth-recovery metrics from the frozen training-derived selected sets;
5. fit the prespecified prediction pipelines using training data;
6. evaluate prediction once on the untouched test cohort.

The proposed test cohort contains 100 samples per phenotype. This may be revised after a runtime and precision pilot, but must be frozen before confirmatory generation.

For real data, an external cohort is preferred. If none exists, repeated nested cross-validation is required: every preprocessing, selection, tuning, and gate operation occurs inside the outer training split, while the outer validation split estimates prediction only.

### 6.3 Methods under comparison

The primary operating-point comparison uses:

- DESeq2 at adjusted `p <= 0.05`;
- edgeR at adjusted `p <= 0.05`;
- Wilcoxon with BH-adjusted `p <= 0.05`;
- the locked, gated SVC-RFECV procedure, where “locked” may include a prespecified training-only tuning rule rather than identical realised hyperparameters in every replicate.

Each method begins from the same training-derived gene universe, while retaining the count representation required by its model.

A secondary gene-budget analysis compares all methods at prespecified values such as `k = 10, 25, 50, 100, 250`. Conventional methods rank by finite unadjusted p-value and then gene identifier; unavailable values appear last. SVC-RFE is fitted to return exactly `k` features.

The gene-budget analysis forces every method to select `k` genes. It assesses ranking only and does not produce supported discoveries. A comparison using the primary SVC set size as `k` must be labelled **performance at the SVC gene budget**, not a symmetric size-matched analysis.

### 6.4 Simulation regimes

The proposed initial regimes are:

| Regime | Truth prevalence | Absolute log2 fold change |
| --- | ---: | ---: |
| Complete null | 0 | 0 |
| Current compcodeR mean-shift | 10% | locked `testset_1` generator distribution |
| Sparse, strong | 1% | 2.0 |
| Dense, weak | 20% | 0.5 |

Truth count is `ceiling(prevalence * number_of_genes)`. Non-null regimes assign approximately half the effects in each direction. Truth genes are sampled across baseline-expression and dispersion strata rather than being associated with leading row positions.

The compcodeR version, every generator argument, and each replicate's latent parameters must be retained. Proposed prevalence and effect values may be changed during pilot work if they create ceiling or floor performance, but must be frozen before confirmation.

The initial grid retains:

- `11 reference + 9 alternative`;
- `20 reference + 16 alternative`;
- `30 reference + 24 alternative`;
- depths `1e5`, `1e6`, and `1e7`.

Run the complete null and current compcodeR regime across that grid. Run sparse/strong and dense/weak at the central `20+16`, depth `1e6` setting. At depth `1e6`, add a fixed-total imbalance comparison of `20+20`, `30+10`, and `10+30`.

Every distinct training design receives a matching complete-null cell with the same sample sizes, depth, dispersion model, filtering path, and test size. Null calibration must not be borrowed from another balance or sample-size design.

Later extensions may consider correlated truth blocks, label noise, technical confounding, cohort shift, and another simulator. Correlated scenarios report gene-level and signal-block recovery separately.

## 7. Benchmark outcomes

### 7.1 Truth-recovery metrics

For selected set `S`, truth set `T`, and universe `U`:

```text
TP = number of genes in both S and T
FP = number of genes in S but not T
FN = number of genes in T but not S
TN = number of genes in neither S nor T

Power = TP / (TP + FN)
FPR = FP / (FP + TN)
FDP = FP / (TP + FP), or 0 when no genes are selected
Youden's index = Power - FPR
```

Power and Youden's index are undefined when the truth set is empty. Youden's index is undefined when there are no background genes. Negative Youden values are retained.

Across independent datasets, report mean power, empirical FDR, mean FPR, mean Youden's index, selected-set size, probability of no selection, and process failure rate. Empirical FDR is the mean replicate-level FDP, not pooled false positives divided by pooled calls.

These metrics can be shown exploratorily in Stage 1, but only Stage 3 supplies the design and replication needed for confirmatory interpretation.

### 7.2 Complete-null behaviour

When there are no truth genes, report:

- probability of selecting any gene;
- expected number of selected genes;
- per-gene false-selection frequency;
- probability that the model-level gate passes;
- independent-test balanced accuracy in Stage 3.

Under the complete null, replicate FDP is one whenever any gene is selected and zero otherwise. Therefore:

```text
empirical FDR = probability of selecting any gene
```

Ungated RFECV must select at least one feature and consequently has null empirical FDR equal to one. This is why candidate features cannot automatically be treated as supported discoveries.

### 7.3 Selection stability

Stage 3 assesses stability at representative scenario cells using 50 stratified bootstrap resamples with replacement within phenotype, preserving the original group sizes. Report:

- gene selection frequencies;
- median pairwise Jaccard similarity among pairs whose union is nonempty;
- proportions of pairs where both sets are empty or exactly one is empty;
- variability in selected-set size;
- truth and background selection frequencies separately.

Do not assign two empty sets Jaccard similarity one and include that value in the median. Report sample-resampling instability separately from CV-partition instability. The latter uses 20 prespecified split seeds on an unchanged dataset.

Stability is descriptive unless a separately specified stability-selection method provides formal guarantees.

### 7.4 Predictive assessment

Independent-test balanced accuracy is the primary predictive metric. Report phenotype sensitivity and specificity alongside it; AUROC is secondary.

To compare the predictive utility of selected gene sets, pass every method's selected genes to the same fixed downstream classifier. The proposal is L2 logistic regression with `C = 1` and balanced class weights after training-derived `log2(CPM + 1)` transformation and scaling. If no genes are selected, use a prespecified no-information classifier rather than dropping the replicate.

The native fitted SVC may also be evaluated, but that answers the separate question of how well the complete SVC procedure predicts. The internal RFECV score remains a tuning statistic and must never be presented as independent predictive performance.

## 8. Estimation and confirmatory decisions

The independently generated dataset—not genes, CV folds, or algorithm seeds—is the unit of uncertainty. Methods are paired within datasets, so method differences retain the replicate pairing.

Report method means, Monte Carlo standard errors, 95% intervals, and paired intervals. Except for binomial event rates, use a seeded nonparametric bootstrap with 10,000 resamples of complete replicate rows. Use a two-sided exact binomial interval for null probability of any selection.

Retryable infrastructure failures rerun the identical dataset and seeds and retain an attempt log. A deterministic analytical failure on valid data remains a failure and must not be replaced, treated as an empty set, or silently omitted. A method with an analytical failure cannot receive a confirmatory compatibility or comparative-performance label for that scenario.

Seven replicates remain suitable for smoke testing and the existing exploratory synthetic plot. The proposed Stage 3 minimum is 200 non-null replicates for the primary cell, increased before confirmation if the pilot indicates that more are needed for a 95% interval half-width of at most `0.02`. The proposed matching-null minimum is 700 replicates. Secondary cells use at least 100 replicates and remain descriptive or inconclusive unless they attain the same precision.

The provisional target empirical FDR is `q = 0.05`, with non-inferiority limit `q + delta_FDP = 0.07`, where `delta_FDP = 0.02`.

A procedure is **empirically false-selection-compatible** only if:

1. the upper 95% confidence bound for null probability of any selection is below `0.07`;
2. the upper 95% confidence bound for non-null empirical FDR is below `0.07`;
3. it has no deterministic analytical failures in that scenario or matching null;
4. its procedure was frozen before confirmatory results were examined.

This is non-inferiority to 7%, not evidence that FDR is below 5%. The stronger claim requires the upper confidence bound below `0.05`.

The proposed definition of **useful power** is mean power at least `0.20` with a lower 95% bound above `0.10`. A procedure is **empirically DEG-compatible** only when it is both false-selection-compatible and usefully powered.

A procedure is **predictive** only when the lower 95% bound for independent-test balanced accuracy exceeds `0.5`.

For paired SVC-minus-comparator differences, the proposed margins are:

- power non-inferiority: lower bound greater than `-0.05`;
- FDP non-inferiority: upper bound below `0.02`;
- material power gain: lower bound greater than `0.05`;
- material FDP reduction: upper bound below `-0.02`;
- stability non-inferiority at common gene budget: lower Jaccard bound greater than `-0.05`.

The proposed primary confirmatory contrast is SVC versus DESeq2 in the compcodeR mean-shift regime at `20+16` training samples and depth `1e6`, with its matching null. Other methods and cells are secondary unless simultaneous intervals or multiplicity correction are prespecified.

All numerical tolerances in this section are proposals for scientific agreement before Stage 3, not settings that need to be settled before Stage 1.

## 9. Stage 3 outputs and regression coverage

### 9.1 Required outputs

The expanded benchmark retains:

- per-replicate truth identities and latent simulation parameters;
- complete candidate and supported gene sets;
- gate statistics, permutation scores, and p-values;
- preprocessing, model, CV, and permutation parameters and seeds;
- training and test group sizes;
- held-out predictions and predictive metrics;
- stability frequencies where calculated;
- process and retry statuses;
- per-replicate metrics, not only summaries;
- the complete fitted prediction objects.

Outputs from the expanded workflow go to a distinct benchmark directory and must not overwrite ordinary pipeline or existing full-synthetic results.

### 9.2 Additional regression coverage

Stage 2 tests verify that:

- a clearly null fixture can return `no_predictive_signal`;
- every gate permutation repeats feature-count optimisation;
- gate outputs distinguish candidate and supported features;
- the existing constrained-permutation results are unchanged;
- a held-out dataset, where supplied for testing, cannot influence the gate.

Stage 3 tests verify that:

- training and test samples are disjoint and correctly aligned;
- filtering, transformations, selection, tuning, and the gate cannot use test samples;
- every requested group size and imbalance direction is preserved;
- every signal scenario has its required matching-null design;
- the end-to-end and conditional truth universes are calculated separately;
- conventional methods and SVC receive the intended shared training universe;
- per-replicate outputs, failure statuses, summaries, and retained plots agree.

## 10. Implementation and review order

1. Agree the Stage 1 reference defaults and the validation rules for parameter overrides and future tuning.
2. Implement Stage 1 changes only in the SVC path and its output readers.
3. Add the opt-in SVC smoke test and confirm that the default non-SVC smoke remains unchanged.
4. Run the rehabilitated SVC against the existing full-synthetic benchmark.
5. Review its candidate-feature behaviour and runtime, then decide whether a separately tested tunable mode is needed.
6. Implement the Stage 2 model-level gate and targeted null checks as a separate change.
7. Review supported-set behaviour and decide whether Stage 3 is justified.
8. If justified, implement the expanded benchmark as an opt-in workflow with separate parameters and outputs.
9. Freeze the complete procedure after development pilots, then generate fresh confirmatory simulations.

Enabling or restoring `mlstep` must not implicitly enable Stage 3. No claim that SVC is a valid DEG alternative follows from Stages 1 or 2 alone. Conversely, Stage 3 should not be implemented merely to complete the mechanical SVC fixes.

## 11. References

- Guyon I, Weston J, Barnhill S, Vapnik V. [Gene Selection for Cancer Classification using Support Vector Machines](https://doi.org/10.1023/A:1012487302797). *Machine Learning*. 2002;46:389–422.
- Ambroise C, McLachlan GJ. [Selection bias in gene extraction on the basis of microarray gene-expression data](https://pubmed.ncbi.nlm.nih.gov/11983868/). *PNAS*. 2002;99:6562–6566.
- Meinshausen N, Bühlmann P. [Stability selection](https://stat.ethz.ch/Manuscripts/buhlmann/stabilityselection.pdf). *Journal of the Royal Statistical Society: Series B*. 2010;72:417–473.
- scikit-learn. [RFECV documentation](https://scikit-learn.org/stable/modules/generated/sklearn.feature_selection.RFECV.html).
- scikit-learn. [Nested versus non-nested cross-validation](https://scikit-learn.org/stable/auto_examples/model_selection/plot_nested_cross_validation_iris.html).
