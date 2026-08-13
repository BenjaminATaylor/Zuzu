#!/usr/bin/env python3

"""Validate the artifacts produced by the opt-in SVC smoke run."""

import json
import os
import pickle
import sys
from argparse import Namespace
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

sys.path.insert(0, "/zuzu-bin")
from zuzu_svc import (
    build_coarse_grid,
    build_refinement_grid,
    choose_smallest_max,
    fit_selected_route,
    fit_svc,
    parse_coarse_grid,
    predict_raw_counts,
)


PRODUCTION_COARSE_COUNTS = [
    1, 2, 3, 5, 8, 13, 20, 32, 50, 80, 125,
    200, 320, 500, 800, 1250, 2000, 3200, 5000, 8000, 10000,
]
SMOKE_COARSE_GRID = "1,10,100"
SMOKE_REFINEMENT_POINTS = 4


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def read_integer_column(frame, column, prefix):
    values = pd.to_numeric(frame[column], errors="raise")
    check(np.isfinite(values).all(), f"{prefix}: {column} contains non-finite values")
    check(np.equal(values, np.floor(values)).all(), f"{prefix}: {column} is not integral")
    return values.astype(int)


def read_pipe_integers(value, prefix):
    if pd.isna(value):
        return []
    try:
        return [int(item) for item in str(value).split("|") if item != ""]
    except ValueError as error:
        raise AssertionError(f"{prefix}: invalid integer schedule") from error


def read_pipe_floats(value, prefix):
    if pd.isna(value):
        return []
    try:
        values = [float(item) for item in str(value).split("|") if item != ""]
    except ValueError as error:
        raise AssertionError(f"{prefix}: invalid score schedule") from error
    check(np.isfinite(values).all(), f"{prefix}: score schedule is non-finite")
    return values


def validate_gate_null_searches(gate_null, coarse_anchors,
                                refinement_points, n_features, prefix):
    required_columns = {
        "coarse_feature_counts",
        "refined_feature_counts",
        "evaluated_feature_counts",
        "mean_test_scores",
        "score_sums",
        "coarse_winner",
        "refinement_lower",
        "refinement_upper",
        "selected_feature_count",
        "selected_stage",
    }
    check(required_columns.issubset(gate_null.columns),
          f"{prefix}: gate null adaptive-search fields are incomplete")
    expected_coarse = list(build_coarse_grid(n_features, coarse_anchors))
    for row_number, row in gate_null.iterrows():
        row_prefix = f"{prefix} permutation {row_number + 1}"
        coarse = read_pipe_integers(row["coarse_feature_counts"], row_prefix)
        refined = read_pipe_integers(row["refined_feature_counts"], row_prefix)
        evaluated = read_pipe_integers(row["evaluated_feature_counts"], row_prefix)
        scores = read_pipe_floats(row["mean_test_scores"], row_prefix)
        score_sums = read_pipe_floats(row["score_sums"], row_prefix)
        winner = int(row["coarse_winner"])
        selected = int(row["selected_feature_count"])
        check(coarse == expected_coarse,
              f"{row_prefix}: coarse search differs from the configured rule")
        expected_refined, expected_lower, expected_upper = build_refinement_grid(
            coarse, winner, refinement_points
        )
        check(refined == list(expected_refined),
              f"{row_prefix}: refinement differs from its own coarse winner")
        check(evaluated == sorted(coarse + refined)
              and len(evaluated) == len(scores) == len(score_sums),
              f"{row_prefix}: evaluated counts and scores are misaligned")
        expected_winner = choose_smallest_max(
            coarse, [score_sums[evaluated.index(count)] for count in coarse]
        )
        check(winner == expected_winner,
              f"{row_prefix}: coarse winner is not its coarse-score maximum")
        check(np.isclose(max(scores), float(row["max_mean_cv_score"])),
              f"{row_prefix}: retained null maximum is incorrect")
        check(selected == choose_smallest_max(evaluated, score_sums),
              f"{row_prefix}: selected count is not its own adaptive maximum")
        expected_stage = "refined" if selected in refined else "coarse"
        check(row["selected_stage"] == expected_stage,
              f"{row_prefix}: selected search stage is incorrect")
        check(int(row["refinement_lower"]) == expected_lower
              and int(row["refinement_upper"]) == expected_upper,
              f"{row_prefix}: refinement bounds are incorrect")


# Check the production schedule without fitting a 10,000-gene model
production_anchors = parse_coarse_grid(
    ",".join(str(count) for count in PRODUCTION_COARSE_COUNTS)
)
check(list(production_anchors) == PRODUCTION_COARSE_COUNTS,
      "production coarse-grid text was not parsed exactly")
check(list(build_coarse_grid(10000, production_anchors))
      == PRODUCTION_COARSE_COUNTS,
      "10,000-feature production coarse grid changed")
production_refined, production_lower, production_upper = build_refinement_grid(
    PRODUCTION_COARSE_COUNTS, 200, 8
)
check(list(production_refined) == [139, 154, 171, 190, 211, 234, 260, 288]
      and production_lower == 125 and production_upper == 320,
      "production refinement grid around 200 features changed")
check(choose_smallest_max([1, 2, 3], [0.5, 1.0, 1.0]) == 2,
      "adaptive search does not resolve score ties toward fewer features")


# Check that a refined selection is reproduced on the complete dataset
route_generator = np.random.default_rng(418)
route_counts = route_generator.normal(size=(20, 100))
route_labels = np.asarray([0] * 10 + [1] * 10)
route_counts[route_labels == 1, :5] += 1.5
route_arguments = Namespace(svc_c=1.0)
route_estimator, route_support, route_ranking, selected_route = fit_selected_route(
    route_counts,
    route_labels,
    coarse_counts=[1, 10, 100],
    refined_counts=[3, 6, 16, 40],
    refinement_upper=100,
    selected_count=3,
    selected_stage="refined",
    arguments=route_arguments,
    resolved_class_weight="balanced",
)
check(selected_route == [100, 40, 16, 6, 3],
      "the complete-data refined route is incorrect")
check(int(route_support.sum()) == 3
      and route_estimator.n_features_in_ == 3,
      "the refined route did not fit the requested feature count")
check(np.equal(route_ranking[route_support], 1).all()
      and np.greater(route_ranking[~route_support], 1).all(),
      "the refined route support and rankings disagree")


# Verify one retained adaptive-search curve against its configured rule
def validate_adaptive_curve(curve, n_features, coarse_anchors,
                            refinement_points, selected_count, prefix):
    required_columns = {
        "n_features",
        "mean_test_score",
        "std_test_score",
        "search_stage",
        "coarse_winner",
        "refinement_lower",
        "refinement_upper",
        "selected_feature_count",
    }
    check(required_columns.issubset(curve.columns),
          f"{prefix}: adaptive-search curve columns are incomplete")

    feature_counts = read_integer_column(curve, "n_features", prefix).tolist()
    check(feature_counts == sorted(set(feature_counts)),
          f"{prefix}: evaluated feature counts are not unique and increasing")
    expected_coarse = list(build_coarse_grid(n_features, coarse_anchors))
    coarse_rows = curve[curve["search_stage"] == "coarse"]
    refined_rows = curve[curve["search_stage"] == "refined"]
    check(set(curve["search_stage"]) <= {"coarse", "refined"},
          f"{prefix}: curve contains an unknown search stage")
    check(read_integer_column(
        coarse_rows, "n_features", f"{prefix} coarse"
    ).tolist() == expected_coarse,
          f"{prefix}: coarse feature-count schedule is incorrect")

    split_score_columns = [
        column for column in curve if column.startswith("split")
        and column.endswith("_test_score")
    ]
    check(len(split_score_columns) == 5,
          f"{prefix}: five fold scores were not retained")
    split_scores = curve[split_score_columns].apply(pd.to_numeric, errors="raise")
    means = pd.to_numeric(curve["mean_test_score"], errors="raise")
    standard_deviations = pd.to_numeric(
        curve["std_test_score"], errors="raise"
    )
    check(np.isfinite(split_scores.to_numpy()).all()
          and np.isfinite(means).all()
          and means.between(0, 1).all(),
          f"{prefix}: balanced-accuracy scores are invalid")
    check(np.allclose(means, split_scores.mean(axis=1))
          and np.allclose(standard_deviations, split_scores.std(axis=1, ddof=0)),
          f"{prefix}: retained score summaries disagree with fold scores")

    score_sums = np.sum(split_scores.to_numpy(dtype=float), axis=1)
    coarse_score_sums = score_sums[coarse_rows.index].tolist()
    expected_coarse_winner = choose_smallest_max(
        expected_coarse, coarse_score_sums
    )
    expected_refined, expected_lower, expected_upper = build_refinement_grid(
        expected_coarse, expected_coarse_winner, refinement_points
    )
    expected_refined = list(expected_refined)
    check(read_integer_column(
        refined_rows, "n_features", f"{prefix} refined"
    ).tolist() == expected_refined,
          f"{prefix}: refined feature-count schedule is incorrect")
    check(feature_counts == sorted(expected_coarse + expected_refined),
          f"{prefix}: curve does not contain the complete adaptive search")

    for column, expected in (
        ("coarse_winner", expected_coarse_winner),
        ("refinement_lower", expected_lower),
        ("refinement_upper", expected_upper),
    ):
        values = read_integer_column(curve, column, prefix)
        check(values.eq(expected).all(),
              f"{prefix}: recorded {column} is incorrect")

    expected_selected = choose_smallest_max(
        feature_counts, score_sums
    )
    chosen = read_integer_column(curve, "selected_feature_count", prefix)
    check(chosen.sum() == 1
          and feature_counts[chosen.tolist().index(1)] == expected_selected
          and selected_count == expected_selected,
          f"{prefix}: final adaptive-search selection is incorrect")
    selected_stage = "refined" if expected_selected in expected_refined else "coarse"
    return {
        "coarse_feature_counts": expected_coarse,
        "refined_feature_counts": expected_refined,
        "evaluated_feature_counts": feature_counts,
        "selected_feature_count": expected_selected,
        "selected_stage": selected_stage,
        "coarse_winner": expected_coarse_winner,
        "refinement_lower": expected_lower,
        "refinement_upper": expected_upper,
        "mean_test_scores": means.tolist(),
        "score_sums": score_sums.tolist(),
    }


# Load the raw and filtered matrices used by the smoke workflow
validation_root = Path("/validation")
artifact_dir = validation_root / "artifacts"
sheet_dir = validation_root / "sheets"
n_permutations = int(sys.argv[1])
reference_label = sys.argv[2]
alternative_label = sys.argv[3]

raw_counts = pd.read_csv(validation_root / "raw_counts.tsv", sep="\t", index_col=0)
analysis_counts = pd.read_csv(validation_root / "analysis_counts.csv", index_col=0)
raw_counts.index = raw_counts.index.astype(str)
raw_counts.columns = raw_counts.columns.astype(str)
analysis_counts.index = analysis_counts.index.astype(str)
analysis_counts.columns = analysis_counts.columns.astype(str)
analysis_genes = analysis_counts.index.tolist()
expected_library_sizes = np.rint(raw_counts).sum(axis=0)
check(len(analysis_genes) > 100, "pare-mode fixture does not contain more than 100 genes")
smoke_coarse_anchors = parse_coarse_grid(SMOKE_COARSE_GRID)

prefixes = [("svc_basic", "basic")]
prefixes.extend(
    (f"svc_permute_replicate_{replicate}", "permutation")
    for replicate in range(1, n_permutations + 1)
)


# Check every baseline and constrained-permutation SVC artifact set
for prefix, expected_analysis_type in prefixes:
    table = pd.read_csv(artifact_dir / f"{prefix}_table.csv")
    curve = pd.read_csv(artifact_dir / f"{prefix}_curve.csv")
    folds = pd.read_csv(artifact_dir / f"{prefix}_folds.csv")
    metadata = json.loads((artifact_dir / f"{prefix}_metadata.json").read_text())
    selected_count = int(
        (artifact_dir / f"{prefix}_n_selected_features.txt").read_text().strip()
    )
    supported_count_text = (
        artifact_dir / f"{prefix}_n_supported_features.txt"
    ).read_text().strip()
    gate = pd.read_csv(artifact_dir / f"{prefix}_gate.csv")
    gate_null = pd.read_csv(artifact_dir / f"{prefix}_gate_null.csv")

    expected_columns = [
        "gene",
        "modelled_feature",
        "selected_feature",
        "supported_feature",
        "selection_rank",
        "coefficient",
        "coefficient_sign",
    ]
    check(table.columns.tolist() == expected_columns, f"{prefix}: unexpected table columns")
    table["gene"] = table["gene"].astype(str)
    check(table["gene"].tolist() == analysis_genes, f"{prefix}: gene universe changed")
    modelled = read_integer_column(table, "modelled_feature", prefix)
    selected = read_integer_column(table, "selected_feature", prefix)
    check(set(modelled.unique()).issubset({0, 1}), f"{prefix}: modelled status is not binary")
    check(set(selected.unique()).issubset({0, 1}), f"{prefix}: selected status is not binary")
    check(modelled.tolist() == [1] * 100 + [0] * (len(table) - 100),
          f"{prefix}: pare mode did not preserve the complete universe")
    check((selected <= modelled).all() and selected.sum() == selected_count,
          f"{prefix}: selected-feature counts disagree")

    # The baseline is gated, while constrained diagnostic permutations are not
    if expected_analysis_type == "basic":
        supported_count = int(supported_count_text)
        supported = read_integer_column(table, "supported_feature", prefix)
        check((supported <= selected).all() and supported.sum() == supported_count,
              f"{prefix}: supported-feature counts disagree")
        check(gate.loc[0, "status"] == "passed",
              f"{prefix}: strong smoke signal did not pass the gate")
        check(supported.tolist() == selected.tolist(),
              f"{prefix}: passing gate did not endorse every candidate")
    else:
        supported_count = None
        check(supported_count_text == "NA",
              f"{prefix}: diagnostic permutation reported a supported count")
        check(table["supported_feature"].isna().all(),
              f"{prefix}: diagnostic permutation endorsed candidate features")
        check(gate.loc[0, "status"] == "not_run_diagnostic_permutation",
              f"{prefix}: constrained permutation unexpectedly ran the gate")
        check(gate_null.empty,
              f"{prefix}: constrained permutation retained gate null runs")

    search_expected = validate_adaptive_curve(
        curve,
        n_features=100,
        coarse_anchors=smoke_coarse_anchors,
        refinement_points=SMOKE_REFINEMENT_POINTS,
        selected_count=selected_count,
        prefix=prefix,
    )
    feature_counts = search_expected["evaluated_feature_counts"]
    means = pd.to_numeric(curve["mean_test_score"], errors="raise")

    sample_sheet = pd.read_csv(sheet_dir / f"{prefix}.csv")
    sample_sheet["sample"] = sample_sheet["sample"].astype(str)
    sample_sheet["phenotype"] = sample_sheet["phenotype"].astype(str)
    expected_mapping = {reference_label: 0, alternative_label: 1}
    expected_encoded = sample_sheet["phenotype"].map(expected_mapping).astype(int)
    check(folds["sample"].astype(str).tolist() == sample_sheet["sample"].tolist(),
          f"{prefix}: fold allocation changed sample order")
    check(read_integer_column(folds, "encoded_class", prefix).tolist()
          == expected_encoded.tolist(), f"{prefix}: class encoding is incorrect")
    check(set(read_integer_column(folds, "validation_fold", prefix)) == {1, 2, 3, 4, 5},
          f"{prefix}: five validation folds were not retained")

    parameters = metadata["parameters"]
    expected_parameters = {
        "kernel": "linear",
        "C": 1.0,
        "class_weight": "balanced",
        "scoring": "balanced_accuracy",
        "folds": 5,
        "shuffle": True,
        "seed": 1729,
        "coarse_grid": list(smoke_coarse_anchors),
        "refinement_points": SMOKE_REFINEMENT_POINTS,
        "min_features_to_select": 1,
        "n_jobs": 1,
    }
    check(metadata["analysis_type"] == expected_analysis_type,
          f"{prefix}: analysis type is incorrect")
    check(metadata["class_mapping"] == expected_mapping,
          f"{prefix}: metadata class mapping is incorrect")
    check(metadata["group_counts"] == {reference_label: 30, alternative_label: 24},
          f"{prefix}: imbalanced group sizes changed")
    check(metadata["library_size_source"]
          == "CLEANINPUTS pre-filter rounded count matrix",
          f"{prefix}: library sizes did not come from the pre-filter matrix")
    check(parameters == expected_parameters, f"{prefix}: resolved parameters are incorrect")
    adaptive_search = metadata["adaptive_search"]
    for key in (
        "coarse_feature_counts",
        "refined_feature_counts",
        "evaluated_feature_counts",
        "selected_feature_count",
        "selected_stage",
        "coarse_winner",
        "refinement_lower",
        "refinement_upper",
    ):
        check(adaptive_search[key] == search_expected[key],
              f"{prefix}: adaptive-search metadata {key} is incorrect")
    check(np.allclose(
        adaptive_search["mean_test_scores"], search_expected["mean_test_scores"]
    ) and np.allclose(
        adaptive_search["score_sums"], search_expected["score_sums"]
    ), f"{prefix}: adaptive-search metadata scores are incorrect")
    check(metadata["estimator_steps"] == ["scaler", "svc"],
          f"{prefix}: estimator steps are incorrect")
    check(metadata["feature_counts"] == {
        "normalization": len(raw_counts),
        "analysis": len(analysis_counts),
        "modelled": 100,
        "selected": selected_count,
        "supported": supported_count,
    }, f"{prefix}: feature-universe counts are incorrect")
    recorded_sizes = pd.Series(metadata["library_sizes"], dtype=float)
    check(recorded_sizes.index.tolist() == sample_sheet["sample"].tolist(),
          f"{prefix}: library-size sample order changed")
    check(np.allclose(recorded_sizes, expected_library_sizes.loc[recorded_sizes.index]),
          f"{prefix}: pre-filter library sizes are incorrect")

    # Check that scaling remains inside every fitted model and score raw counts
    with (artifact_dir / f"{prefix}_model.pickle").open("rb") as handle:
        bundle = pickle.load(handle)
    check(bundle["normalization_gene_ids"] == raw_counts.index.tolist(),
          f"{prefix}: normalization universe was not saved")
    check(bundle["analysis_gene_ids"] == analysis_genes,
          f"{prefix}: analysis universe was not saved")
    check(bundle["selected_gene_ids"] == table.loc[selected == 1, "gene"].tolist(),
          f"{prefix}: selected-gene order was not saved")
    expected_supported_genes = None if supported_count is None else table.loc[
        table["supported_feature"] == 1, "gene"
    ].tolist()
    check(bundle["supported_gene_ids"] == expected_supported_genes,
          f"{prefix}: supported-gene order was not saved")
    check(bundle["signal_gate_passed"] == (expected_analysis_type == "basic"),
          f"{prefix}: saved-model support status is incorrect")
    check(bundle["schema_version"] == "zuzu_svc_model_v3",
          f"{prefix}: saved model does not use the adaptive-search schema")
    check(bundle["search"] == adaptive_search,
          f"{prefix}: saved-model adaptive search disagrees with metadata")
    check(isinstance(bundle["estimator"], Pipeline),
          f"{prefix}: fitted estimator is not a pipeline")
    check(isinstance(bundle["estimator"].named_steps["scaler"], StandardScaler),
          f"{prefix}: scaler is outside the fitted estimator")
    check(isinstance(bundle["estimator"].named_steps["svc"], SVC),
          f"{prefix}: fitted classifier is not an SVC")
    predictions = np.asarray(predict_raw_counts(bundle, raw_counts))
    check(predictions.ndim == 1 and len(predictions) == len(sample_sheet),
          f"{prefix}: saved bundle could not score every raw-count sample")
    check(set(predictions).issubset({0, 1}),
          f"{prefix}: saved bundle returned an unknown class")

    # Independently reconstruct training log-CPM and compare saved-model predictions
    rounded_counts = raw_counts.apply(pd.to_numeric, errors="raise").round()
    training_library_sizes = rounded_counts.loc[
        bundle["normalization_gene_ids"]
    ].sum(axis=0)
    model_log_cpm = np.log2(
        rounded_counts.loc[bundle["model_gene_ids"]].divide(
            training_library_sizes, axis=1
        ) * 1e6 + 1.0
    ).transpose()
    expected_predictions = np.asarray(
        bundle["estimator"].predict(
            model_log_cpm.loc[:, bundle["selected_gene_ids"]].to_numpy(dtype=float)
        )
    )
    check(np.array_equal(predictions, expected_predictions),
          f"{prefix}: saved raw-count preprocessing changed fitted predictions")

    if expected_analysis_type == "basic":
        # Check the exact add-one gate calculation and every adaptive null search
        check(len(gate) == 1 and len(gate_null) == 19,
              f"{prefix}: gate did not retain all requested null searches")
        observed_score = float(means.max())
        null_scores = pd.to_numeric(gate_null["max_mean_cv_score"], errors="raise")
        null_counts = read_integer_column(
            gate_null, "selected_feature_count", f"{prefix} gate"
        )
        exceedances = int((null_scores >= observed_score).sum())
        check(np.isclose(float(gate.loc[0, "observed_score"]), observed_score),
              f"{prefix}: gate statistic is not the adaptive-search maximum")
        check((null_counts >= 1).all() and (null_counts <= 100).all(),
              f"{prefix}: gate permutation selected an invalid feature count")
        validate_gate_null_searches(
            gate_null,
            coarse_anchors=smoke_coarse_anchors,
            refinement_points=SMOKE_REFINEMENT_POINTS,
            n_features=100,
            prefix=f"{prefix} gate",
        )
        check(read_integer_column(
            gate_null, "at_least_observed", f"{prefix} gate"
        ).tolist() == (null_scores >= observed_score).astype(int).tolist(),
              f"{prefix}: gate tie/exceedance rule is incorrect")
        check(int(gate.loc[0, "exceedances"]) == exceedances,
              f"{prefix}: gate exceedance count is incorrect")
        check(np.isclose(float(gate.loc[0, "p_value"]), (1 + exceedances) / 20),
              f"{prefix}: gate add-one p-value is incorrect")
        check(gate_null["permuted_encoded_classes"].nunique() == len(gate_null),
              f"{prefix}: gate label allocations are not unique")
        observed_classes = "|".join(folds["encoded_class"].astype(str).tolist())
        null_classes = gate_null["permuted_encoded_classes"].astype(str)
        null_allocations = null_classes.map(lambda value: value.split("|"))
        check((null_classes != observed_classes).all(),
              f"{prefix}: gate reused the observed label allocation")
        check(null_allocations.map(len).eq(len(folds)).all()
              and null_allocations.map(lambda value: set(value) <= {"0", "1"}).all()
              and null_allocations.map(lambda value: value.count("1")).eq(
                  int(folds["encoded_class"].sum())
              ).all(),
              f"{prefix}: gate permutations did not preserve phenotype margins")


# Check the opposite imbalance direction and the unsupported small-class error
with TemporaryDirectory(dir="/tmp") as focused_directory:
    previous_directory = os.getcwd()
    os.chdir(focused_directory)
    try:
        arguments = Namespace(
            samplesheet=str(sheet_dir / "svc_basic.csv"),
            counts=str(validation_root / "analysis_counts.csv"),
            library_sizes=None,
            library_genes=None,
            reference=alternative_label,
            pare="true",
            svc_c=1.0,
            class_weight="balanced",
            folds=5,
            scoring="balanced_accuracy",
            seed=1729,
            coarse_grid=parse_coarse_grid(SMOKE_COARSE_GRID),
            refinement_points=SMOKE_REFINEMENT_POINTS,
            n_jobs=1,
            gate="false",
            gate_permutations=19,
            gate_alpha=0.05,
            gate_seed=2718,
            gate_early_stop="true",
            output_prefix="svc_reference_minority",
            analysis_type="focused_reference_minority",
            replicate=None,
            depth=None,
            reference_samples=None,
            alternative_samples=None,
        )
        fit_svc(arguments)
        minority_metadata = json.loads(
            Path("svc_reference_minority_metadata.json").read_text()
        )
        check(minority_metadata["class_mapping"]
              == {alternative_label: 0, reference_label: 1},
              "reference-minority class direction is incorrect")
        check(minority_metadata["group_counts"]
              == {alternative_label: 24, reference_label: 30},
              "reference-minority group sizes changed")

        basic_sheet = pd.read_csv(sheet_dir / "svc_basic.csv")
        small_sheet = pd.concat([
            basic_sheet[basic_sheet["phenotype"] == reference_label].head(4),
            basic_sheet[basic_sheet["phenotype"] == alternative_label].head(6),
        ], ignore_index=True)
        small_sheet.to_csv("small_sheet.csv", index=False)
        analysis_counts.loc[:, small_sheet["sample"]].to_csv("small_counts.csv")
        small_arguments = Namespace(**{
            **vars(arguments),
            "samplesheet": "small_sheet.csv",
            "counts": "small_counts.csv",
            "reference": reference_label,
            "output_prefix": "svc_small_class",
        })
        try:
            fit_svc(small_arguments)
        except ValueError as error:
            check("requires at least 5 samples in each phenotype" in str(error)
                  and f"{reference_label}=4" in str(error),
                  "small-class error did not explain the unsupported design")
        else:
            raise AssertionError("SVC accepted a phenotype group smaller than its fold count")

        # A deterministic complete-null fixture must retain candidates only as diagnostics
        null_samples = [f"null_sample_{index}" for index in range(1, 11)]
        null_sheet = pd.DataFrame(
            {
                "sample": null_samples,
                "phenotype": [reference_label] * 5 + [alternative_label] * 5,
            }
        )
        null_counts = pd.DataFrame(
            np.repeat(np.arange(100, 112, dtype=int)[:, None], 10, axis=1),
            index=[f"null_gene_{index}" for index in range(1, 13)],
            columns=null_samples,
        )
        null_sheet.to_csv("null_sheet.csv", index=False)
        null_counts.to_csv("null_counts.csv")
        null_arguments = Namespace(**{
            **vars(arguments),
            "samplesheet": "null_sheet.csv",
            "counts": "null_counts.csv",
            "reference": reference_label,
            "pare": "false",
            "gate": "true",
            "output_prefix": "svc_complete_null_1",
            "analysis_type": "focused_complete_null",
        })
        fit_svc(null_arguments)
        fit_svc(Namespace(**{
            **vars(null_arguments),
            "output_prefix": "svc_complete_null_2",
        }))

        null_table = pd.read_csv("svc_complete_null_1_table.csv")
        null_gate = pd.read_csv("svc_complete_null_1_gate.csv")
        null_gate_scores = pd.read_csv("svc_complete_null_1_gate_null.csv")
        null_metadata = json.loads(
            Path("svc_complete_null_1_metadata.json").read_text()
        )
        check(null_gate.loc[0, "status"] == "no_predictive_signal",
              "complete-null fixture passed the SVC signal gate")
        check(bool(null_gate.loc[0, "early_stopped"])
              and int(null_gate.loc[0, "permutations_run"]) == 1,
              "complete-null gate did not stop after failure became certain")
        check(pd.isna(null_gate.loc[0, "p_value"])
              and null_gate.loc[0, "p_value_statement"] == "p > 0.05",
              "early-stopped gate reported an invalid exact p-value")
        check(int(null_gate.loc[0, "exceedances"]) == 1
              and np.isclose(null_gate.loc[0, "p_value_lower_bound"], 0.1),
              "complete-null gate lower bound is incorrect")
        check(len(null_gate_scores) == 1
              and int(null_gate_scores.loc[0, "at_least_observed"]) == 1,
              "complete-null gate did not count an equal null score conservatively")
        validate_gate_null_searches(
            null_gate_scores,
            coarse_anchors=smoke_coarse_anchors,
            refinement_points=SMOKE_REFINEMENT_POINTS,
            n_features=12,
            prefix="complete-null gate",
        )
        check(int(null_table["selected_feature"].sum()) >= 1
              and int(null_table["supported_feature"].sum()) == 0,
              "failed gate did not separate forced candidates from supported features")
        check(Path("svc_complete_null_1_n_supported_features.txt").read_text().strip()
              == "0", "failed gate did not publish a zero supported count")
        check(null_metadata["supported_gene_ids"] == []
              and not null_metadata["signal_gate"]["status"] == "passed",
              "complete-null metadata incorrectly endorses a signature")

        # The separate gate seed must reproduce both null allocations and decisions
        check(Path("svc_complete_null_1_gate.csv").read_text()
              == Path("svc_complete_null_2_gate.csv").read_text(),
              "gate summary is not reproducible from its recorded seeds")
        check(Path("svc_complete_null_1_gate_null.csv").read_text()
              == Path("svc_complete_null_2_gate_null.csv").read_text(),
              "gate null allocation is not reproducible from its recorded seeds")
        with Path("svc_complete_null_1_model.pickle").open("rb") as handle:
            null_bundle = pickle.load(handle)
        try:
            predict_raw_counts(null_bundle, null_counts, require_supported=True)
        except ValueError as error:
            check("signal gate did not pass" in str(error),
                  "unsupported prediction error is unclear")
        else:
            raise AssertionError("failed-gate model was accepted as supported prediction")
    finally:
        os.chdir(previous_directory)

print("SVC smoke artifact, gate, class-direction, and raw-prediction checks passed.")
