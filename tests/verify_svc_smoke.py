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
from zuzu_svc import fit_svc, predict_raw_counts


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def read_integer_column(frame, column, prefix):
    values = pd.to_numeric(frame[column], errors="raise")
    check(np.isfinite(values).all(), f"{prefix}: {column} contains non-finite values")
    check(np.equal(values, np.floor(values)).all(), f"{prefix}: {column} is not integral")
    return values.astype(int)


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

    expected_columns = [
        "gene",
        "modelled_feature",
        "selected_feature",
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

    feature_counts = read_integer_column(curve, "n_features", prefix).tolist()
    means = pd.to_numeric(curve["mean_test_score"], errors="raise")
    chosen = read_integer_column(curve, "selected_feature_count", prefix)
    split_score_columns = [
        column for column in curve if column.startswith("split")
        and column.endswith("_test_score")
    ]
    check(feature_counts == [1, 50, 100], f"{prefix}: actual RFECV counts are wrong")
    check(len(split_score_columns) == 5, f"{prefix}: five fold scores were not retained")
    check(np.isfinite(means).all() and means.between(0, 1).all(),
          f"{prefix}: balanced-accuracy scores are invalid")
    check(chosen.sum() == 1 and feature_counts[chosen.tolist().index(1)] == selected_count,
          f"{prefix}: selected curve point disagrees with the table")
    smallest_maximum = min(
        count for count, score in zip(feature_counts, means) if score == means.max()
    )
    check(selected_count == smallest_maximum,
          f"{prefix}: score ties did not select the smallest feature set")

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
        "rfe_step": 50,
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
    check(metadata["estimator_steps"] == ["scaler", "svc"],
          f"{prefix}: estimator steps are incorrect")
    check(metadata["feature_counts"] == {
        "normalization": len(raw_counts),
        "analysis": len(analysis_counts),
        "modelled": 100,
        "selected": selected_count,
    }, f"{prefix}: feature-universe counts are incorrect")
    recorded_sizes = pd.Series(metadata["library_sizes"], dtype=float)
    check(recorded_sizes.index.tolist() == sample_sheet["sample"].tolist(),
          f"{prefix}: library-size sample order changed")
    check(np.allclose(recorded_sizes, expected_library_sizes.loc[recorded_sizes.index]),
          f"{prefix}: pre-filter library sizes are incorrect")

    # Check that scaling remains inside RFECV and the saved bundle predicts raw counts
    with (artifact_dir / f"{prefix}_model.pickle").open("rb") as handle:
        bundle = pickle.load(handle)
    check(bundle["normalization_gene_ids"] == raw_counts.index.tolist(),
          f"{prefix}: normalization universe was not saved")
    check(bundle["analysis_gene_ids"] == analysis_genes,
          f"{prefix}: analysis universe was not saved")
    check(bundle["selected_gene_ids"] == table.loc[selected == 1, "gene"].tolist(),
          f"{prefix}: selected-gene order was not saved")
    check(isinstance(bundle["rfecv"].estimator, Pipeline),
          f"{prefix}: RFECV estimator is not a pipeline")
    check(isinstance(bundle["rfecv"].estimator.named_steps["scaler"], StandardScaler),
          f"{prefix}: scaler is outside RFECV")
    check(isinstance(bundle["rfecv"].estimator.named_steps["svc"], SVC),
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
        bundle["rfecv"].predict(model_log_cpm)
    )
    check(np.array_equal(predictions, expected_predictions),
          f"{prefix}: saved raw-count preprocessing changed fitted predictions")


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
            step=50,
            n_jobs=1,
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
    finally:
        os.chdir(previous_directory)

print("SVC smoke artifact, class-direction, and raw-prediction checks passed.")
