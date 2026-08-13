#!/usr/bin/env python3

"""Fit and inspect Zuzu's shared adaptive SVC RFE-CV procedure."""

import argparse
import json
import math
import os
import pickle
import platform
import re
import tempfile
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.metrics import get_scorer
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC


ALLOWED_SCORING = {
    "accuracy",
    "balanced_accuracy",
    "f1",
    "f1_macro",
    "f1_weighted",
    "roc_auc",
}

DEFAULT_COARSE_GRID = (
    1,
    2,
    3,
    5,
    8,
    13,
    20,
    32,
    50,
    80,
    125,
    200,
    320,
    500,
    800,
    1250,
    2000,
    3200,
    5000,
    8000,
    10000,
)
DEFAULT_COARSE_GRID_TEXT = ",".join(str(value) for value in DEFAULT_COARSE_GRID)


# Read and validate the two aligned analysis inputs
def read_analysis_inputs(samplesheet_path, counts_path, reference_label):
    samplesheet = pd.read_csv(samplesheet_path)
    required_columns = {"sample", "phenotype"}
    missing_columns = required_columns.difference(samplesheet.columns)
    if missing_columns:
        raise ValueError(
            "Sample sheet is missing columns: " + ", ".join(sorted(missing_columns))
        )

    if samplesheet[["sample", "phenotype"]].isna().any().any():
        raise ValueError("Sample sheet contains missing sample or phenotype values")
    samplesheet = samplesheet.copy()
    samplesheet["sample"] = samplesheet["sample"].astype(str)
    samplesheet["phenotype"] = samplesheet["phenotype"].astype(str)
    if (samplesheet["sample"].str.strip() == "").any():
        raise ValueError("Sample sheet contains blank sample identifiers")
    if (samplesheet["phenotype"].str.strip() == "").any():
        raise ValueError("Sample sheet contains blank phenotype values")
    if samplesheet["sample"].duplicated().any():
        raise ValueError("Sample sheet contains duplicate sample identifiers")

    counts = pd.read_csv(counts_path, index_col=0)
    if counts.index.isna().any():
        raise ValueError("Count matrix contains missing gene identifiers")
    counts.index = counts.index.astype(str)
    counts.columns = counts.columns.astype(str)
    if counts.index.duplicated().any():
        raise ValueError("Count matrix contains duplicate gene identifiers")
    if counts.columns.duplicated().any():
        raise ValueError("Count matrix contains duplicate sample identifiers")

    sample_order = samplesheet["sample"].tolist()
    if counts.columns.tolist() != sample_order:
        raise ValueError("Count columns do not match sample sheet order")
    if counts.shape[0] < 2:
        raise ValueError("SVC RFE-CV requires at least two genes")

    try:
        counts = counts.apply(pd.to_numeric, errors="raise")
    except (TypeError, ValueError) as error:
        raise ValueError("Count matrix contains non-numeric values") from error
    count_values = counts.to_numpy(dtype=float)
    if not np.isfinite(count_values).all():
        raise ValueError("Count matrix contains missing or non-finite values")
    if (count_values < 0).any():
        raise ValueError("Count matrix contains negative values")

    phenotype_levels = samplesheet["phenotype"].unique().tolist()
    if len(phenotype_levels) != 2:
        raise ValueError("SVC requires exactly two phenotypes")
    if reference_label not in phenotype_levels:
        raise ValueError(
            f"Requested reference level '{reference_label}' is absent from the sample sheet"
        )
    alternative_levels = [
        level for level in phenotype_levels if level != reference_label
    ]
    alternative_label = alternative_levels[0]

    encoded = np.where(
        samplesheet["phenotype"].to_numpy() == reference_label,
        0,
        1,
    )
    return samplesheet, counts, encoded, alternative_label


# Read the pre-filter library sizes supplied by CLEANINPUTS
def read_normalization_sidecars(
    library_sizes_path,
    library_genes_path,
    sample_order,
    analysis_gene_ids,
):
    library_sizes_frame = pd.read_csv(library_sizes_path)
    if library_sizes_frame.columns.tolist() != ["sample", "library_size"]:
        raise ValueError(
            "SVC library-size sidecar must contain sample and library_size columns"
        )
    if library_sizes_frame.isna().any().any():
        raise ValueError("SVC library-size sidecar contains missing values")
    library_sizes_frame["sample"] = library_sizes_frame["sample"].astype(str)
    if library_sizes_frame["sample"].duplicated().any():
        raise ValueError("SVC library-size sidecar contains duplicate samples")
    if library_sizes_frame["sample"].tolist() != sample_order:
        raise ValueError("SVC library sizes do not match sample sheet order")
    try:
        library_sizes = pd.to_numeric(
            library_sizes_frame["library_size"], errors="raise"
        ).to_numpy(dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError("SVC library sizes must be numeric") from error
    if not np.isfinite(library_sizes).all() or (library_sizes <= 0).any():
        raise ValueError("SVC library sizes must be positive and finite")

    library_genes_frame = pd.read_csv(library_genes_path)
    if library_genes_frame.columns.tolist() != ["gene"]:
        raise ValueError("SVC library-gene sidecar must contain one gene column")
    if library_genes_frame["gene"].isna().any():
        raise ValueError("SVC library-gene sidecar contains missing identifiers")
    library_gene_ids = library_genes_frame["gene"].astype(str).tolist()
    if len(library_gene_ids) != len(set(library_gene_ids)):
        raise ValueError("SVC library-gene sidecar contains duplicate identifiers")
    missing_analysis_genes = set(analysis_gene_ids).difference(library_gene_ids)
    if missing_analysis_genes:
        raise ValueError(
            "SVC analysis genes are absent from the normalization universe: "
            + ", ".join(sorted(missing_analysis_genes)[:5])
        )

    return library_sizes, library_gene_ids


# Calculate library sizes directly for the existing full-synthetic input
def derive_normalization(counts):
    library_sizes = counts.sum(axis=0).to_numpy(dtype=float)
    if not np.isfinite(library_sizes).all() or (library_sizes <= 0).any():
        raise ValueError("SVC normalization requires every library size to be positive")
    return library_sizes, counts.index.tolist()


# Convert the filtered analysis counts to log2 counts per million
def calculate_log_cpm(counts, library_sizes):
    scaled_counts = counts.divide(library_sizes, axis="columns") * 1e6
    transformed = np.log2(scaled_counts + 1.0).transpose()
    if not np.isfinite(transformed.to_numpy(dtype=float)).all():
        raise ValueError("SVC log2(CPM + 1) transformation produced non-finite values")
    return transformed


# Parse the tunable reference grid supplied through Nextflow
def parse_coarse_grid(value):
    try:
        feature_counts = [int(item.strip()) for item in value.split(",")]
    except (AttributeError, TypeError, ValueError) as error:
        raise argparse.ArgumentTypeError(
            "--coarse-grid must be a comma-separated list of integers"
        ) from error
    if len(feature_counts) < 2:
        raise argparse.ArgumentTypeError(
            "--coarse-grid must contain at least two feature counts"
        )
    if feature_counts[0] != 1:
        raise argparse.ArgumentTypeError("--coarse-grid must start at 1")
    if any(value < 1 for value in feature_counts):
        raise argparse.ArgumentTypeError(
            "--coarse-grid feature counts must all be positive"
        )
    if any(left >= right for left, right in zip(feature_counts, feature_counts[1:])):
        raise argparse.ArgumentTypeError(
            "--coarse-grid feature counts must be strictly increasing"
        )
    return feature_counts


# Map the reference grid to the available feature universe
def build_coarse_grid(input_feature_count, reference_grid):
    input_feature_count = int(input_feature_count)
    reference_grid = np.asarray(reference_grid, dtype=int)
    if input_feature_count < 2:
        raise ValueError("Adaptive SVC RFE-CV requires at least two features")
    if len(reference_grid) < 2 or reference_grid[0] != 1:
        raise ValueError("The coarse reference grid must start at 1")
    if np.any(np.diff(reference_grid) <= 0):
        raise ValueError("The coarse reference grid must be strictly increasing")

    # Search very small feature universes exhaustively
    if input_feature_count <= len(reference_grid):
        return np.arange(1, input_feature_count + 1, dtype=int)

    reference_total = int(reference_grid[-1])
    if reference_total <= 1:
        raise ValueError("The coarse reference grid must end above 1")
    relative_log_positions = np.log(reference_grid) / math.log(reference_total)
    mapped = np.rint(
        np.exp(math.log(input_feature_count) * relative_log_positions)
    ).astype(int)
    mapped = np.clip(mapped, 1, input_feature_count)
    feature_counts = np.unique(
        np.concatenate((mapped, np.asarray([1, input_feature_count], dtype=int)))
    )
    return feature_counts.astype(int)


# Resolve score ties in favour of the smallest evaluated feature set
def choose_smallest_max(feature_counts, score_sums):
    feature_counts = np.asarray(feature_counts, dtype=int)
    score_sums = np.asarray(score_sums, dtype=float)
    if feature_counts.shape != score_sums.shape or feature_counts.ndim != 1:
        raise ValueError("Feature counts and score sums must be aligned vectors")
    if len(feature_counts) == 0 or not np.isfinite(score_sums).all():
        raise ValueError("Feature-count selection requires finite scores")
    maximum = np.max(score_sums)
    return int(np.min(feature_counts[score_sums == maximum]))


# Generate new integer counts around one run's own coarse maximum
def build_refinement_grid(coarse_counts, coarse_winner, refinement_points):
    coarse_counts = np.asarray(coarse_counts, dtype=int)
    coarse_winner = int(coarse_winner)
    refinement_points = int(refinement_points)
    winner_positions = np.flatnonzero(coarse_counts == coarse_winner)
    if len(winner_positions) != 1:
        raise ValueError("The coarse winner must occur exactly once in the grid")
    winner_position = int(winner_positions[0])
    lower = int(
        coarse_counts[winner_position - 1]
        if winner_position > 0
        else coarse_winner
    )
    upper = int(
        coarse_counts[winner_position + 1]
        if winner_position < len(coarse_counts) - 1
        else coarse_winner
    )
    if refinement_points == 0 or lower == upper:
        return np.asarray([], dtype=int), lower, upper

    proposed = np.rint(
        np.geomspace(lower, upper, refinement_points + 2)[1:-1]
    ).astype(int)
    coarse_set = set(int(value) for value in coarse_counts)
    refined = sorted(
        {
            int(value)
            for value in proposed
            if lower < value < upper and int(value) not in coarse_set
        }
    )
    return np.asarray(refined, dtype=int), lower, upper


# Record which validation fold contains each sample
def make_fold_table(cv, transformed_counts, encoded, samplesheet):
    validation_fold = np.zeros(len(samplesheet), dtype=int)
    for fold_number, (_, validation_indices) in enumerate(
        cv.split(transformed_counts, encoded), start=1
    ):
        validation_fold[validation_indices] = fold_number
    if (validation_fold == 0).any():
        raise RuntimeError("Could not assign every sample to a validation fold")

    return pd.DataFrame(
        {
            "sample": samplesheet["sample"],
            "phenotype": samplesheet["phenotype"],
            "encoded_class": encoded,
            "validation_fold": validation_fold,
        }
    )


# Convert values used by numpy and pandas into JSON-safe Python values
def json_safe(value):
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return [json_safe(item) for item in value.tolist()]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, Path):
        return str(value)
    return value


# Validate a raw count frame before applying a saved model bundle
def prepare_raw_counts(bundle, raw_counts):
    if not isinstance(raw_counts, pd.DataFrame):
        raise TypeError("Raw counts must be supplied as a pandas DataFrame")
    if raw_counts.index.isna().any():
        raise ValueError("Raw count matrix contains missing gene identifiers")

    raw_counts = raw_counts.copy()
    raw_counts.index = raw_counts.index.astype(str)
    raw_counts.columns = raw_counts.columns.astype(str)
    if raw_counts.index.duplicated().any():
        raise ValueError("Raw count matrix contains duplicate gene identifiers")
    if raw_counts.columns.duplicated().any():
        raise ValueError("Raw count matrix contains duplicate sample identifiers")

    normalization_gene_ids = bundle["normalization_gene_ids"]
    selected_gene_ids = bundle["selected_gene_ids"]
    missing_genes = set(normalization_gene_ids).difference(raw_counts.index)
    if missing_genes:
        raise ValueError(
            "Raw count matrix is missing genes required for normalization: "
            + ", ".join(sorted(missing_genes)[:5])
        )

    try:
        numeric_counts = raw_counts.loc[normalization_gene_ids].apply(
            pd.to_numeric, errors="raise"
        )
    except (TypeError, ValueError) as error:
        raise ValueError("Raw count matrix contains non-numeric values") from error
    numeric_values = numeric_counts.to_numpy(dtype=float)
    if not np.isfinite(numeric_values).all():
        raise ValueError("Raw count matrix contains missing or non-finite values")
    if (numeric_values < 0).any():
        raise ValueError("Raw count matrix contains negative values")

    # Match CLEANINPUTS, which rounds counts before calculating library sizes
    rounded_counts = pd.DataFrame(
        np.rint(numeric_values),
        index=numeric_counts.index,
        columns=numeric_counts.columns,
    )
    library_sizes = rounded_counts.sum(axis=0)
    if (library_sizes <= 0).any():
        raise ValueError("Raw count matrix contains an empty sample library")

    selected_counts = rounded_counts.loc[selected_gene_ids]
    transformed = np.log2(selected_counts.divide(library_sizes, axis=1) * 1e6 + 1.0)
    return transformed.transpose()


# Score new raw-count samples using all saved preprocessing information
def predict_raw_counts(bundle, raw_counts, require_supported=False):
    if require_supported and not bundle.get("signal_gate_passed"):
        raise ValueError("The SVC signal gate did not pass for this signature")
    transformed = prepare_raw_counts(bundle, raw_counts)
    return bundle["estimator"].predict(transformed.to_numpy(dtype=float))


# Return the fitted SVC decision function for new raw-count samples
def decision_function_raw_counts(bundle, raw_counts, require_supported=False):
    if require_supported and not bundle.get("signal_gate_passed"):
        raise ValueError("The SVC signal gate did not pass for this signature")
    transformed = prepare_raw_counts(bundle, raw_counts)
    return bundle["estimator"].decision_function(
        transformed.to_numpy(dtype=float)
    )


# Construct a fresh scaler and linear classifier for one elimination fit
def make_estimator(arguments, resolved_class_weight):
    return Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "svc",
                SVC(
                    kernel="linear",
                    C=arguments.svc_c,
                    class_weight=resolved_class_weight,
                ),
            ),
        ]
    )


# Rank active features by squared coefficient with stable index tie handling
def get_importances(estimator):
    coefficients = np.asarray(
        estimator.named_steps["svc"].coef_,
        dtype=float,
    )
    importances = np.sum(np.square(coefficients), axis=0)
    if not np.isfinite(importances).all():
        raise RuntimeError("SVC produced missing or non-finite feature importances")
    return importances


# Remove the weakest active features to reach the next scheduled count
def eliminate_to_count(active_features, importances, target_count):
    active_features = np.asarray(active_features, dtype=int)
    importances = np.asarray(importances, dtype=float)
    target_count = int(target_count)
    if len(active_features) != len(importances):
        raise RuntimeError("Could not align active features with SVC importances")
    if target_count < 1 or target_count >= len(active_features):
        raise RuntimeError("The elimination target must be below the active count")

    elimination_order = np.lexsort((active_features, importances))
    remove_count = len(active_features) - target_count
    removed = active_features[elimination_order[:remove_count]]
    retained = active_features[elimination_order[remove_count:]]
    return np.sort(retained), np.sort(removed)


# Score the descending coarse path within one training and validation fold
def score_coarse_fold(
    fold_number,
    training_counts,
    training_labels,
    validation_counts,
    validation_labels,
    descending_counts,
    arguments,
    resolved_class_weight,
    verbose,
):
    active_features = np.arange(training_counts.shape[1], dtype=int)
    scorer = get_scorer(arguments.scoring)
    scores = {}
    states = {}

    for position, feature_count in enumerate(descending_counts):
        feature_count = int(feature_count)
        if len(active_features) != feature_count:
            raise RuntimeError("The coarse RFE path reached an unexpected feature count")
        estimator = make_estimator(arguments, resolved_class_weight)
        estimator.fit(training_counts[:, active_features], training_labels)
        score = float(
            scorer(
                estimator,
                validation_counts[:, active_features],
                validation_labels,
            )
        )
        if not np.isfinite(score):
            raise RuntimeError("Adaptive RFE-CV produced a non-finite fold score")
        importances = get_importances(estimator)
        scores[feature_count] = score
        states[feature_count] = (active_features.copy(), importances.copy())
        if verbose:
            print(
                f"Adaptive RFE-CV coarse fold {fold_number}: "
                f"scored {feature_count} features"
            )

        if position < len(descending_counts) - 1:
            next_count = int(descending_counts[position + 1])
            active_features, _ = eliminate_to_count(
                active_features,
                importances,
                next_count,
            )

    return {"scores": scores, "states": states}


# Continue one fold from its saved upper-neighbour coarse state
def score_refinement_fold(
    fold_number,
    training_counts,
    training_labels,
    validation_counts,
    validation_labels,
    upper_state,
    descending_refined_counts,
    arguments,
    resolved_class_weight,
    verbose,
):
    active_features, importances = upper_state
    active_features = active_features.copy()
    importances = importances.copy()
    scorer = get_scorer(arguments.scoring)
    scores = {}

    for feature_count in descending_refined_counts:
        feature_count = int(feature_count)
        active_features, _ = eliminate_to_count(
            active_features,
            importances,
            feature_count,
        )
        estimator = make_estimator(arguments, resolved_class_weight)
        estimator.fit(training_counts[:, active_features], training_labels)
        score = float(
            scorer(
                estimator,
                validation_counts[:, active_features],
                validation_labels,
            )
        )
        if not np.isfinite(score):
            raise RuntimeError("Adaptive RFE-CV produced a non-finite fold score")
        importances = get_importances(estimator)
        scores[feature_count] = score
        if verbose:
            print(
                f"Adaptive RFE-CV refined fold {fold_number}: "
                f"scored {feature_count} features"
            )

    return scores


# Reproduce the selected coarse or refined route on the complete dataset
def fit_selected_route(
    model_counts,
    encoded,
    coarse_counts,
    refined_counts,
    refinement_upper,
    selected_count,
    selected_stage,
    arguments,
    resolved_class_weight,
):
    descending_coarse = sorted((int(value) for value in coarse_counts), reverse=True)
    if selected_stage == "coarse":
        route = [value for value in descending_coarse if value >= selected_count]
    elif selected_stage == "refined":
        route = [value for value in descending_coarse if value >= refinement_upper]
        route.extend(
            value
            for value in sorted(
                (int(item) for item in refined_counts),
                reverse=True,
            )
            if value >= selected_count
        )
    else:
        raise RuntimeError("The selected adaptive-search stage is invalid")
    if not route or route[-1] != selected_count:
        raise RuntimeError("Could not reconstruct the selected adaptive RFE route")

    model_values = np.asarray(model_counts, dtype=float)
    active_features = np.arange(model_values.shape[1], dtype=int)
    eliminated_batches = []
    selected_estimator = None
    for position, feature_count in enumerate(route):
        if len(active_features) != feature_count:
            raise RuntimeError("The full-data RFE path reached an unexpected count")
        estimator = make_estimator(arguments, resolved_class_weight)
        estimator.fit(model_values[:, active_features], encoded)
        if position == len(route) - 1:
            selected_estimator = estimator
            break
        importances = get_importances(estimator)
        active_features, removed = eliminate_to_count(
            active_features,
            importances,
            route[position + 1],
        )
        eliminated_batches.append(removed)

    support = np.zeros(model_values.shape[1], dtype=bool)
    support[active_features] = True
    ranking = np.ones(model_values.shape[1], dtype=int)
    for position, removed in enumerate(eliminated_batches):
        ranking[removed] = len(eliminated_batches) - position + 1
    return selected_estimator, support, ranking, route


# Fit one complete adaptive search for the observed or a permuted label vector
def fit_adaptive_rfe_cv(
    model_counts,
    encoded,
    arguments,
    resolved_class_weight,
    verbose,
    fit_final,
):
    model_values = np.asarray(model_counts, dtype=float)
    coarse_counts = build_coarse_grid(
        model_values.shape[1],
        arguments.coarse_grid,
    )
    descending_coarse = coarse_counts[::-1]
    cv = StratifiedKFold(
        n_splits=arguments.folds,
        shuffle=True,
        random_state=arguments.seed,
    )
    splits = list(cv.split(model_values, encoded))

    # Run the same recursive coarse path independently inside every fold
    coarse_fold_results = joblib.Parallel(
        n_jobs=arguments.n_jobs,
        prefer="threads",
    )(
        joblib.delayed(score_coarse_fold)(
            fold_number,
            model_values[training_indices],
            encoded[training_indices],
            model_values[validation_indices],
            encoded[validation_indices],
            descending_coarse,
            arguments,
            resolved_class_weight,
            verbose,
        )
        for fold_number, (training_indices, validation_indices) in enumerate(
            splits,
            start=1,
        )
    )
    coarse_split_scores = np.asarray(
        [
            [result["scores"][int(count)] for count in coarse_counts]
            for result in coarse_fold_results
        ],
        dtype=float,
    )
    coarse_score_sums = np.sum(coarse_split_scores, axis=0)
    coarse_winner = choose_smallest_max(coarse_counts, coarse_score_sums)
    refined_counts, refinement_lower, refinement_upper = build_refinement_grid(
        coarse_counts,
        coarse_winner,
        arguments.refinement_points,
    )

    # Branch from each fold's saved upper-neighbour state for local refinement
    if len(refined_counts) > 0:
        descending_refined = refined_counts[::-1]
        refined_fold_scores = joblib.Parallel(
            n_jobs=arguments.n_jobs,
            prefer="threads",
        )(
            joblib.delayed(score_refinement_fold)(
                fold_number,
                model_values[training_indices],
                encoded[training_indices],
                model_values[validation_indices],
                encoded[validation_indices],
                coarse_fold_results[fold_number - 1]["states"][refinement_upper],
                descending_refined,
                arguments,
                resolved_class_weight,
                verbose,
            )
            for fold_number, (training_indices, validation_indices) in enumerate(
                splits,
                start=1,
            )
        )
        refined_split_scores = np.asarray(
            [
                [result[int(count)] for count in refined_counts]
                for result in refined_fold_scores
            ],
            dtype=float,
        )
    else:
        refined_split_scores = np.empty((arguments.folds, 0), dtype=float)

    # Combine the two labelled paths and choose the smallest overall maximum
    evaluated_counts = np.sort(np.concatenate((coarse_counts, refined_counts)))
    coarse_positions = {int(value): index for index, value in enumerate(coarse_counts)}
    refined_positions = {
        int(value): index for index, value in enumerate(refined_counts)
    }
    split_scores = np.empty((arguments.folds, len(evaluated_counts)), dtype=float)
    stages = []
    for position, feature_count in enumerate(evaluated_counts):
        feature_count = int(feature_count)
        if feature_count in coarse_positions:
            split_scores[:, position] = coarse_split_scores[
                :,
                coarse_positions[feature_count],
            ]
            stages.append("coarse")
        else:
            split_scores[:, position] = refined_split_scores[
                :,
                refined_positions[feature_count],
            ]
            stages.append("refined")
    mean_scores = np.mean(split_scores, axis=0)
    std_scores = np.std(split_scores, axis=0)
    score_sums = np.sum(split_scores, axis=0)
    if not (
        np.isfinite(mean_scores).all()
        and np.isfinite(std_scores).all()
        and np.isfinite(score_sums).all()
    ):
        raise RuntimeError("Adaptive RFE-CV produced non-finite summary scores")
    selected_count = choose_smallest_max(evaluated_counts, score_sums)
    selected_position = int(np.flatnonzero(evaluated_counts == selected_count)[0])
    selected_stage = stages[selected_position]

    estimator = None
    support = None
    ranking = None
    selected_route = None
    if fit_final:
        estimator, support, ranking, selected_route = fit_selected_route(
            model_values,
            encoded,
            coarse_counts,
            refined_counts,
            refinement_upper,
            selected_count,
            selected_stage,
            arguments,
            resolved_class_weight,
        )

    search = {
        "procedure": "adaptive two-branch SVC RFE-CV",
        "coarse_reference_grid": list(arguments.coarse_grid),
        "coarse_feature_counts": coarse_counts,
        "refined_feature_counts": refined_counts,
        "evaluated_feature_counts": evaluated_counts,
        "search_stages": stages,
        "coarse_winner": coarse_winner,
        "refinement_lower": refinement_lower,
        "refinement_upper": refinement_upper,
        "selected_feature_count": selected_count,
        "selected_stage": selected_stage,
        "selected_route": selected_route,
        "mean_test_scores": mean_scores,
        "std_test_scores": std_scores,
        "score_sums": score_sums,
        "split_test_scores": split_scores,
        "selection_rule": (
            "smallest feature count attaining the maximum sum of fold scores"
        ),
        "refinement_rule": (
            "logarithmic integer counts inside the neighbours of this run's "
            "smallest coarse maximum"
        ),
        "path_interpretation": (
            "Refined points branch from the saved upper-neighbour state and do not "
            "form one uninterrupted path with coarse points below that state."
        ),
        "importance_rule": (
            "sum of squared linear-SVC coefficients; ties broken by original "
            "feature order"
        ),
    }
    return {
        "search": search,
        "cv": cv,
        "estimator": estimator,
        "support": support,
        "ranking": ranking,
    }


# Generate reproducible, unique label allocations for the SVC signal gate
def iter_gate_permutations(encoded, count, seed):
    observed = tuple(int(value) for value in encoded)
    seen = {observed}
    seed_generator = np.random.default_rng(seed)
    accepted = 0

    while accepted < count:
        permutation_seed = int(
            seed_generator.integers(0, 2**32, dtype=np.uint64)
        )
        permuted = np.random.default_rng(permutation_seed).permutation(encoded)
        signature = tuple(int(value) for value in permuted)
        if signature in seen:
            continue
        seen.add(signature)
        accepted += 1
        yield accepted, permutation_seed, permuted


# Serialize one adaptive search vector for the gate audit table
def pipe_values(values, formatter=str):
    return "|".join(formatter(value) for value in values)


# Test whether the optimised adaptive score exceeds its permuted-label null
def run_signal_gate(
    model_counts,
    encoded,
    arguments,
    resolved_class_weight,
    observed_score,
):
    requested = int(arguments.gate_permutations)
    alpha = float(arguments.gate_alpha)
    maximum_passing_exceedances = math.floor(alpha * (requested + 1) - 1 + 1e-12)
    early_stop = arguments.gate_early_stop == "true"
    exceedances = 0
    null_rows = []

    # Repeat the complete feature-count optimisation for every null allocation
    for permutation, permutation_seed, permuted_encoded in iter_gate_permutations(
        encoded,
        requested,
        arguments.gate_seed,
    ):
        null_result = fit_adaptive_rfe_cv(
            model_counts,
            permuted_encoded,
            arguments,
            resolved_class_weight,
            verbose=0,
            fit_final=False,
        )
        null_search = null_result["search"]
        null_score = float(np.max(null_search["mean_test_scores"]))
        at_least_observed = null_score >= observed_score
        exceedances += int(at_least_observed)
        null_rows.append(
            {
                "permutation": permutation,
                "permutation_seed": permutation_seed,
                "max_mean_cv_score": null_score,
                "selected_feature_count": int(
                    null_search["selected_feature_count"]
                ),
                "selected_stage": null_search["selected_stage"],
                "coarse_winner": int(null_search["coarse_winner"]),
                "refinement_lower": int(null_search["refinement_lower"]),
                "refinement_upper": int(null_search["refinement_upper"]),
                "coarse_feature_counts": pipe_values(
                    null_search["coarse_feature_counts"],
                    lambda value: str(int(value)),
                ),
                "refined_feature_counts": pipe_values(
                    null_search["refined_feature_counts"],
                    lambda value: str(int(value)),
                ),
                "evaluated_feature_counts": pipe_values(
                    null_search["evaluated_feature_counts"],
                    lambda value: str(int(value)),
                ),
                "mean_test_scores": pipe_values(
                    null_search["mean_test_scores"],
                    lambda value: format(float(value), ".17g"),
                ),
                "score_sums": pipe_values(
                    null_search["score_sums"],
                    lambda value: format(float(value), ".17g"),
                ),
                "at_least_observed": int(at_least_observed),
                "permuted_encoded_classes": "|".join(
                    str(int(value)) for value in permuted_encoded
                ),
            }
        )

        # Once too many null scores exceed the observation, passing is impossible
        if early_stop and exceedances > maximum_passing_exceedances:
            break

    permutations_run = len(null_rows)
    stopped_early = permutations_run < requested
    p_value_lower_bound = (1 + exceedances) / (1 + requested)
    if stopped_early:
        p_value = None
        p_value_statement = f"p > {alpha:g}"
        status = "no_predictive_signal"
    else:
        p_value = p_value_lower_bound
        p_value_statement = f"p = {p_value:.12g}"
        status = "passed" if p_value <= alpha else "no_predictive_signal"

    gate_result = {
        "status": status,
        "statistic": "maximum mean adaptive RFE-CV score",
        "scoring": arguments.scoring,
        "observed_score": float(observed_score),
        "p_value": p_value,
        "p_value_statement": p_value_statement,
        "p_value_lower_bound": float(p_value_lower_bound),
        "alpha": alpha,
        "permutations_requested": requested,
        "permutations_run": permutations_run,
        "exceedances": exceedances,
        "early_stopped": stopped_early,
        "early_stop_enabled": early_stop,
        "cv_seed": int(arguments.seed),
        "gate_seed": int(arguments.gate_seed),
        "null_sampling": "unique unrestricted fixed-margin label permutations",
        "tie_rule": "null score >= observed score counts as an exceedance",
        "interpretation": (
            "Model-level test of the global phenotype-label exchangeability null; "
            "not independent predictive performance or gene-wise error control."
        ),
    }
    return gate_result, pd.DataFrame(null_rows)


# Check user-adjustable SVC parameters before fitting
def validate_parameters(arguments):
    if not np.isfinite(arguments.svc_c) or arguments.svc_c <= 0:
        raise ValueError("--svc-c must be a positive finite number")
    if arguments.folds < 2:
        raise ValueError("--folds must be an integer of at least 2")
    try:
        parse_coarse_grid(",".join(str(value) for value in arguments.coarse_grid))
    except argparse.ArgumentTypeError as error:
        raise ValueError(str(error)) from error
    if arguments.refinement_points < 0:
        raise ValueError("--refinement-points cannot be negative")
    if arguments.seed < 0 or arguments.seed > (2**32 - 1):
        raise ValueError("--seed must be between 0 and 4294967295")
    if arguments.gate_seed < 0 or arguments.gate_seed > (2**32 - 1):
        raise ValueError("--gate-seed must be between 0 and 4294967295")
    if arguments.n_jobs == 0:
        raise ValueError("--n-jobs cannot be zero")
    if arguments.scoring not in ALLOWED_SCORING:
        raise ValueError(
            "--scoring must be one of: " + ", ".join(sorted(ALLOWED_SCORING))
        )
    if arguments.class_weight not in {"balanced", "none"}:
        raise ValueError("--class-weight must be either balanced or none")
    if arguments.pare not in {"true", "false"}:
        raise ValueError("--pare must be either true or false")
    if arguments.gate not in {"true", "false"}:
        raise ValueError("--gate must be either true or false")
    if arguments.gate_early_stop not in {"true", "false"}:
        raise ValueError("--gate-early-stop must be either true or false")
    if arguments.gate_permutations < 1:
        raise ValueError("--gate-permutations must be a positive integer")
    if not np.isfinite(arguments.gate_alpha) or not 0 < arguments.gate_alpha < 1:
        raise ValueError("--gate-alpha must be greater than 0 and less than 1")
    if arguments.gate == "true":
        if 1 / (arguments.gate_permutations + 1) > arguments.gate_alpha:
            raise ValueError(
                "--gate-permutations is too small for the requested --gate-alpha; "
                "even zero exceedances could not pass"
            )
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", arguments.output_prefix):
        raise ValueError("--output-prefix contains unsupported characters")


# Fit the corrected SVC procedure and retain all supporting outputs
def fit_svc(arguments):
    validate_parameters(arguments)
    samplesheet, counts, encoded, alternative_label = read_analysis_inputs(
        arguments.samplesheet,
        arguments.counts,
        arguments.reference,
    )
    sample_order = samplesheet["sample"].tolist()
    analysis_gene_ids = counts.index.tolist()

    if bool(arguments.library_sizes) != bool(arguments.library_genes):
        raise ValueError(
            "--library-sizes and --library-genes must be supplied together"
        )
    if arguments.library_sizes:
        library_sizes, normalization_gene_ids = read_normalization_sidecars(
            arguments.library_sizes,
            arguments.library_genes,
            sample_order,
            analysis_gene_ids,
        )
        library_size_source = "CLEANINPUTS pre-filter rounded count matrix"
    else:
        library_sizes, normalization_gene_ids = derive_normalization(counts)
        library_size_source = "all genes in the emitted full-synthetic count matrix"

    group_counts = {
        arguments.reference: int(np.sum(encoded == 0)),
        alternative_label: int(np.sum(encoded == 1)),
    }
    if bool(arguments.reference_samples is not None) != bool(
        arguments.alternative_samples is not None
    ):
        raise ValueError(
            "--reference-samples and --alternative-samples must be supplied together"
        )
    if arguments.reference_samples is not None:
        if arguments.reference_samples < 1 or arguments.alternative_samples < 1:
            raise ValueError("Expected phenotype group sizes must both be positive")
        observed_sizes = (
            group_counts[arguments.reference],
            group_counts[alternative_label],
        )
        expected_sizes = (
            arguments.reference_samples,
            arguments.alternative_samples,
        )
        if observed_sizes != expected_sizes:
            raise ValueError(
                "SVC phenotype counts do not match the requested group sizes: "
                f"observed ref={observed_sizes[0]}, alt={observed_sizes[1]}; "
                f"expected ref={expected_sizes[0]}, alt={expected_sizes[1]}"
            )
    too_small = {
        label: count for label, count in group_counts.items() if count < arguments.folds
    }
    if too_small:
        details = ", ".join(f"{label}={count}" for label, count in too_small.items())
        raise ValueError(
            f"SVC {arguments.folds}-fold cross-validation requires at least "
            f"{arguments.folds} samples in each phenotype; too small: {details}. "
            "Lower --svc_folds only if fewer folds are scientifically appropriate."
        )

    # Apply log-CPM before the fold-specific scaling pipeline
    transformed_counts = calculate_log_cpm(counts, library_sizes)
    pare_mode = arguments.pare == "true"
    model_gene_ids = analysis_gene_ids[:100] if pare_mode else analysis_gene_ids
    model_counts = transformed_counts.loc[:, model_gene_ids]
    if arguments.gate == "true":
        attainable_null_allocations = (
            math.comb(len(encoded), int(np.sum(encoded == 1))) - 1
        )
        if arguments.gate_permutations > attainable_null_allocations:
            raise ValueError(
                "--gate-permutations exceeds the number of unique non-observed "
                f"fixed-margin label allocations ({attainable_null_allocations})"
            )

    # Setup the fold-local scaler, classifier and recursive feature selection
    resolved_class_weight = (
        None if arguments.class_weight == "none" else arguments.class_weight
    )

    # Select the smallest feature set at the maximum fold-score sum
    adaptive_result = fit_adaptive_rfe_cv(
        model_counts,
        encoded,
        arguments,
        resolved_class_weight,
        verbose=2,
        fit_final=True,
    )
    adaptive_search = adaptive_result["search"]
    cv = adaptive_result["cv"]
    feature_counts = np.asarray(
        adaptive_search["evaluated_feature_counts"],
        dtype=int,
    )
    mean_scores = np.asarray(adaptive_search["mean_test_scores"], dtype=float)
    std_scores = np.asarray(adaptive_search["std_test_scores"], dtype=float)
    selected_feature_count = int(adaptive_search["selected_feature_count"])
    selected_indices = np.flatnonzero(adaptive_result["support"])
    selected_gene_ids = [model_gene_ids[index] for index in selected_indices]
    coefficients = np.asarray(
        adaptive_result["estimator"].named_steps["svc"].coef_,
        dtype=float,
    ).reshape(-1)
    if len(coefficients) != len(selected_gene_ids):
        raise RuntimeError("Could not align selected genes with fitted SVC coefficients")

    # Gate the forced candidate set using unrestricted phenotype-label permutations
    observed_score = float(np.max(mean_scores))
    if arguments.gate == "true":
        gate_result, gate_null_table = run_signal_gate(
            model_counts,
            encoded,
            arguments,
            resolved_class_weight,
            observed_score,
        )
        gate_passed = gate_result["status"] == "passed"
        supported_gene_ids = selected_gene_ids if gate_passed else []
    else:
        gate_result = {
            "status": "not_run_diagnostic_permutation",
            "statistic": "maximum mean adaptive RFE-CV score",
            "scoring": arguments.scoring,
            "observed_score": observed_score,
            "p_value": None,
            "p_value_statement": "not run",
            "p_value_lower_bound": None,
            "alpha": float(arguments.gate_alpha),
            "permutations_requested": int(arguments.gate_permutations),
            "permutations_run": 0,
            "exceedances": None,
            "early_stopped": False,
            "early_stop_enabled": arguments.gate_early_stop == "true",
            "cv_seed": int(arguments.seed),
            "gate_seed": int(arguments.gate_seed),
            "null_sampling": "not run",
            "tie_rule": "null score >= observed score counts as an exceedance",
            "interpretation": (
                "Gate deliberately omitted for the constrained-permutation diagnostic."
            ),
        }
        gate_null_table = pd.DataFrame(
            columns=[
                "permutation",
                "permutation_seed",
                "max_mean_cv_score",
                "selected_feature_count",
                "selected_stage",
                "coarse_winner",
                "refinement_lower",
                "refinement_upper",
                "coarse_feature_counts",
                "refined_feature_counts",
                "evaluated_feature_counts",
                "mean_test_scores",
                "score_sums",
                "at_least_observed",
                "permuted_encoded_classes",
            ]
        )
        supported_gene_ids = None

    # Write a complete feature table over the analysis gene universe
    model_position = {gene: index for index, gene in enumerate(model_gene_ids)}
    rank_by_gene = {
        gene: int(adaptive_result["ranking"][index])
        for gene, index in model_position.items()
    }
    coefficient_by_gene = dict(zip(selected_gene_ids, coefficients))
    selected_gene_set = set(selected_gene_ids)
    model_gene_set = set(model_gene_ids)
    output_table = pd.DataFrame({"gene": analysis_gene_ids})
    output_table["modelled_feature"] = output_table["gene"].isin(model_gene_set).astype(int)
    output_table["selected_feature"] = output_table["gene"].isin(selected_gene_set).astype(int)
    if supported_gene_ids is None:
        output_table["supported_feature"] = pd.array(
            [pd.NA] * len(output_table), dtype="Int64"
        )
    else:
        supported_gene_set = set(supported_gene_ids)
        output_table["supported_feature"] = output_table["gene"].isin(
            supported_gene_set
        ).astype(int)
    output_table["selection_rank"] = pd.array(
        [rank_by_gene.get(gene, pd.NA) for gene in analysis_gene_ids], dtype="Int64"
    )
    output_table["coefficient"] = [
        coefficient_by_gene.get(gene, np.nan) for gene in analysis_gene_ids
    ]
    output_table["coefficient_sign"] = pd.array(
        [
            int(np.sign(coefficient_by_gene[gene]))
            if gene in coefficient_by_gene
            else pd.NA
            for gene in analysis_gene_ids
        ],
        dtype="Int64",
    )

    # Write both adaptive search paths and every fold score
    curve_table = pd.DataFrame(
        {
            "n_features": feature_counts,
            "mean_test_score": mean_scores,
            "std_test_score": std_scores,
            "search_stage": adaptive_search["search_stages"],
            "coarse_winner": int(adaptive_search["coarse_winner"]),
            "refinement_lower": int(adaptive_search["refinement_lower"]),
            "refinement_upper": int(adaptive_search["refinement_upper"]),
            "selected_feature_count": (
                feature_counts == selected_feature_count
            ).astype(int),
        }
    )
    for fold_position, fold_scores in enumerate(
        adaptive_search["split_test_scores"],
        start=0,
    ):
        curve_table[f"split{fold_position}_test_score"] = np.asarray(
            fold_scores,
            dtype=float,
        )
    fold_table = make_fold_table(cv, model_counts, encoded, samplesheet)

    output_prefix = arguments.output_prefix
    table_path = f"{output_prefix}_table.csv"
    curve_path = f"{output_prefix}_curve.csv"
    folds_path = f"{output_prefix}_folds.csv"
    plot_path = f"{output_prefix}_fs.png"
    model_path = f"{output_prefix}_model.pickle"
    metadata_path = f"{output_prefix}_metadata.json"
    feature_count_path = f"{output_prefix}_n_selected_features.txt"
    supported_count_path = f"{output_prefix}_n_supported_features.txt"
    gate_path = f"{output_prefix}_gate.csv"
    gate_null_path = f"{output_prefix}_gate_null.csv"
    gate_plot_path = f"{output_prefix}_gate.png"
    output_table.to_csv(table_path, index=False)
    curve_table.to_csv(curve_path, index=False)
    fold_table.to_csv(folds_path, index=False)
    Path(feature_count_path).write_text(
        f"{selected_feature_count}\n",
        encoding="utf-8",
    )
    supported_count = None if supported_gene_ids is None else len(supported_gene_ids)
    Path(supported_count_path).write_text(
        "NA\n" if supported_count is None else f"{supported_count}\n",
        encoding="utf-8",
    )

    # Retain the gate summary and every completed null adaptive search
    gate_summary_table = pd.DataFrame(
        [
            {
                "status": gate_result["status"],
                "statistic": gate_result["statistic"],
                "scoring": gate_result["scoring"],
                "observed_score": gate_result["observed_score"],
                "p_value": gate_result["p_value"],
                "p_value_statement": gate_result["p_value_statement"],
                "p_value_lower_bound": gate_result["p_value_lower_bound"],
                "alpha": gate_result["alpha"],
                "permutations_requested": gate_result["permutations_requested"],
                "permutations_run": gate_result["permutations_run"],
                "exceedances": gate_result["exceedances"],
                "early_stopped": gate_result["early_stopped"],
                "early_stop_enabled": gate_result["early_stop_enabled"],
                "cv_seed": gate_result["cv_seed"],
                "gate_seed": gate_result["gate_seed"],
                "null_sampling": gate_result["null_sampling"],
                "tie_rule": gate_result["tie_rule"],
            }
        ]
    )
    gate_summary_table.to_csv(gate_path, index=False)
    gate_null_table.to_csv(gate_null_path, index=False)

    # Plot the feature-selection curve against real evaluated feature counts
    os.environ.setdefault(
        "MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "zuzu-matplotlib")
    )
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots()
    stage_styles = {
        "coarse": {"marker": "o", "label": "Coarse search"},
        "refined": {"marker": "s", "label": "Refinement branch"},
    }
    search_stages = np.asarray(adaptive_search["search_stages"])
    for stage, style in stage_styles.items():
        stage_mask = search_stages == stage
        if not stage_mask.any():
            continue
        axis.errorbar(
            feature_counts[stage_mask],
            mean_scores[stage_mask],
            yerr=std_scores[stage_mask],
            marker=style["marker"],
            capsize=3,
            linestyle="none",
            label=style["label"],
        )
    selected_score = mean_scores[feature_counts == selected_feature_count][0]
    axis.scatter(
        [selected_feature_count],
        [selected_score],
        color="red",
        zorder=3,
        label="Selected feature count",
    )
    if feature_counts[-1] / feature_counts[0] >= 10:
        axis.set_xscale("log")
    axis.set_xlabel("Number of SVC candidate features")
    axis.set_ylabel(f"Mean cross-validation score ({arguments.scoring})")
    axis.legend()
    figure.tight_layout()
    figure.savefig(plot_path, format="png", dpi=150)
    plt.close(figure)

    # Plot the gate null distribution separately from the feature-selection curve
    gate_figure, gate_axis = plt.subplots()
    if len(gate_null_table) > 0:
        gate_axis.hist(
            gate_null_table["max_mean_cv_score"],
            bins=min(20, max(5, len(gate_null_table))),
            color="white",
            edgecolor="black",
        )
    gate_axis.axvline(
        observed_score,
        color="red",
        label="Observed adaptive RFE-CV maximum",
    )
    gate_axis.set_xlabel(
        f"Maximum mean adaptive RFE-CV score ({arguments.scoring})"
    )
    gate_axis.set_ylabel("Number of label permutations")
    gate_axis.set_title(f"Signal gate: {gate_result['status']}")
    gate_progress = (
        f"{gate_result['permutations_run']}/"
        f"{gate_result['permutations_requested']} permutations"
    )
    if gate_result["early_stopped"]:
        gate_progress += "\nEarly-stopped partial null"
    gate_axis.text(
        0.02,
        0.98,
        gate_progress,
        transform=gate_axis.transAxes,
        horizontalalignment="left",
        verticalalignment="top",
    )
    gate_axis.legend()
    gate_figure.tight_layout()
    gate_figure.savefig(gate_plot_path, format="png", dpi=150)
    plt.close(gate_figure)

    resolved_parameters = {
        "kernel": "linear",
        "C": float(arguments.svc_c),
        "class_weight": resolved_class_weight,
        "scoring": arguments.scoring,
        "folds": int(arguments.folds),
        "shuffle": True,
        "seed": int(arguments.seed),
        "coarse_grid": list(arguments.coarse_grid),
        "refinement_points": int(arguments.refinement_points),
        "min_features_to_select": 1,
        "n_jobs": int(arguments.n_jobs),
    }
    preprocessing = {
        "input": "raw counts",
        "training_library_sizes": library_size_source,
        "new_sample_library_sizes": "stored normalization gene universe",
        "round_new_counts": True,
        "normalization": "counts per million",
        "transform": "log2(CPM + 1)",
        "cpm_scale": 1e6,
        "pseudocount": 1.0,
        "scaling": "StandardScaler fitted inside each RFE-CV training fold",
    }
    class_mapping = {arguments.reference: 0, alternative_label: 1}
    versions = {
        "python": platform.python_version(),
        "scikit_learn": sklearn.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "pandas": pd.__version__,
        "joblib": joblib.__version__,
    }

    # Save a model bundle that can be loaded without custom pickle classes
    serialized_search = json_safe(adaptive_search)
    model_bundle = {
        "schema_version": "zuzu_svc_model_v3",
        "search": serialized_search,
        "estimator": adaptive_result["estimator"],
        "support": adaptive_result["support"],
        "ranking": adaptive_result["ranking"],
        "normalization_gene_ids": normalization_gene_ids,
        "analysis_gene_ids": analysis_gene_ids,
        "model_gene_ids": model_gene_ids,
        "selected_gene_ids": selected_gene_ids,
        "supported_gene_ids": supported_gene_ids,
        "signal_gate_passed": gate_result["status"] == "passed",
        "gate": gate_result,
        "class_mapping": class_mapping,
        "preprocessing": preprocessing,
        "parameters": resolved_parameters,
        "software_versions": versions,
    }
    with open(model_path, "wb") as model_handle:
        pickle.dump(model_bundle, model_handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Save human-readable run metadata alongside the fitted model
    metadata = {
        "schema_version": "zuzu_svc_metadata_v3",
        "analysis_type": arguments.analysis_type,
        "context": {
            "replicate": arguments.replicate,
            "depth": arguments.depth,
            "reference_samples": arguments.reference_samples,
            "alternative_samples": arguments.alternative_samples,
        },
        "sample_order": sample_order,
        "class_mapping": class_mapping,
        "group_counts": group_counts,
        "library_sizes": dict(zip(sample_order, library_sizes)),
        "library_size_source": library_size_source,
        "parameters": resolved_parameters,
        "preprocessing": preprocessing,
        "feature_counts": {
            "normalization": len(normalization_gene_ids),
            "analysis": len(analysis_gene_ids),
            "modelled": len(model_gene_ids),
            "selected": len(selected_gene_ids),
            "supported": supported_count,
        },
        "selected_gene_ids": selected_gene_ids,
        "supported_gene_ids": supported_gene_ids,
        "adaptive_search": serialized_search,
        "cv_feature_counts": feature_counts,
        "cv_mean_test_scores": mean_scores,
        "cv_std_test_scores": std_scores,
        "estimator_steps": list(adaptive_result["estimator"].named_steps),
        "pare_mode": pare_mode,
        "signal_gate": gate_result,
        "coefficient_interpretation": (
            "Signs describe the fitted multivariable SVC boundary for alternative "
            "class 1 versus reference class 0; they are not expression-effect directions."
        ),
        "files": {
            "table": table_path,
            "curve": curve_path,
            "folds": folds_path,
            "plot": plot_path,
            "model": model_path,
            "selected_feature_count": feature_count_path,
            "supported_feature_count": supported_count_path,
            "gate": gate_path,
            "gate_null": gate_null_path,
            "gate_plot": gate_plot_path,
        },
        "software_versions": versions,
    }
    with open(metadata_path, "w", encoding="utf-8") as metadata_handle:
        json.dump(json_safe(metadata), metadata_handle, indent=2, allow_nan=False)
        metadata_handle.write("\n")

    print(
        f"Selected {selected_feature_count} of {len(model_gene_ids)} "
        "modelled features"
    )
    if supported_count is None:
        print("Signal gate not run for the constrained-permutation diagnostic")
    else:
        print(
            f"Signal gate {gate_result['status']}; supported {supported_count} features"
        )


# Define the command-line interface shared by all three Nextflow modules
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--library-sizes")
    parser.add_argument("--library-genes")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--pare", required=True)
    parser.add_argument("--svc-c", required=True, type=float)
    parser.add_argument("--class-weight", required=True)
    parser.add_argument("--folds", required=True, type=int)
    parser.add_argument("--scoring", required=True)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument(
        "--coarse-grid",
        required=True,
        type=parse_coarse_grid,
    )
    parser.add_argument("--refinement-points", required=True, type=int)
    parser.add_argument("--n-jobs", required=True, type=int)
    parser.add_argument("--gate", required=True)
    parser.add_argument("--gate-permutations", required=True, type=int)
    parser.add_argument("--gate-alpha", required=True, type=float)
    parser.add_argument("--gate-seed", required=True, type=int)
    parser.add_argument("--gate-early-stop", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--analysis-type", required=True)
    parser.add_argument("--replicate", type=int)
    parser.add_argument("--depth", type=float)
    parser.add_argument("--reference-samples", type=int)
    parser.add_argument("--alternative-samples", type=int)
    return parser.parse_args()


def main():
    fit_svc(parse_arguments())


if __name__ == "__main__":
    main()
