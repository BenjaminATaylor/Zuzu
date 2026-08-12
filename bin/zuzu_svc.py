#!/usr/bin/env python3

"""Fit and inspect Zuzu's shared SVC-RFECV procedure."""

import argparse
import json
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
from sklearn.feature_selection import RFECV
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
        raise ValueError("SVC-RFECV requires at least two genes")

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


# Reconstruct evaluated feature counts for older scikit-learn releases
def get_feature_counts(rfecv, input_feature_count, step, minimum):
    if "n_features" in rfecv.cv_results_:
        feature_counts = np.asarray(rfecv.cv_results_["n_features"], dtype=int)
    else:
        descending_counts = [input_feature_count]
        while descending_counts[-1] > minimum:
            descending_counts.append(max(minimum, descending_counts[-1] - step))
        feature_counts = np.asarray(descending_counts[::-1], dtype=int)

    score_count = len(rfecv.cv_results_["mean_test_score"])
    if len(feature_counts) != score_count:
        raise RuntimeError(
            "Could not align RFECV scores with their evaluated feature counts"
        )
    if feature_counts[0] != minimum or feature_counts[-1] != input_feature_count:
        raise RuntimeError("RFECV returned an unexpected feature-count sequence")
    if len(feature_counts) > 1 and not np.all(np.diff(feature_counts) > 0):
        raise RuntimeError("RFECV feature counts are not strictly increasing")
    return feature_counts


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
def predict_raw_counts(bundle, raw_counts):
    transformed = prepare_raw_counts(bundle, raw_counts)
    return bundle["estimator"].predict(transformed.to_numpy(dtype=float))


# Return the fitted SVC decision function for new raw-count samples
def decision_function_raw_counts(bundle, raw_counts):
    transformed = prepare_raw_counts(bundle, raw_counts)
    return bundle["estimator"].decision_function(
        transformed.to_numpy(dtype=float)
    )


# Check user-adjustable SVC parameters before fitting
def validate_parameters(arguments):
    if not np.isfinite(arguments.svc_c) or arguments.svc_c <= 0:
        raise ValueError("--svc-c must be a positive finite number")
    if arguments.folds < 2:
        raise ValueError("--folds must be an integer of at least 2")
    if arguments.step < 1:
        raise ValueError("--step must be a positive integer")
    if arguments.seed < 0 or arguments.seed > (2**32 - 1):
        raise ValueError("--seed must be between 0 and 4294967295")
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

    # Setup the fold-local scaler, classifier and recursive feature selection
    resolved_class_weight = (
        None if arguments.class_weight == "none" else arguments.class_weight
    )
    estimator = Pipeline(
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
    cv = StratifiedKFold(
        n_splits=arguments.folds,
        shuffle=True,
        random_state=arguments.seed,
    )
    rfecv = RFECV(
        estimator=estimator,
        importance_getter="named_steps.svc.coef_",
        step=arguments.step,
        cv=cv,
        scoring=arguments.scoring,
        min_features_to_select=1,
        n_jobs=arguments.n_jobs,
        verbose=2,
    )

    # Select the smallest feature set at the maximum mean fold score
    rfecv.fit(model_counts, encoded)
    feature_counts = get_feature_counts(
        rfecv,
        input_feature_count=len(model_gene_ids),
        step=arguments.step,
        minimum=1,
    )
    mean_scores = np.asarray(rfecv.cv_results_["mean_test_score"], dtype=float)
    std_scores = np.asarray(rfecv.cv_results_["std_test_score"], dtype=float)
    if not np.isfinite(mean_scores).all() or not np.isfinite(std_scores).all():
        raise RuntimeError("RFECV produced missing or non-finite validation scores")
    tied_best_indices = np.flatnonzero(mean_scores == np.max(mean_scores))
    expected_selected_count = int(feature_counts[tied_best_indices[0]])
    if int(rfecv.n_features_) != expected_selected_count:
        raise RuntimeError(
            "RFECV did not choose the smallest feature set at the maximum score"
        )

    selected_indices = rfecv.get_support(indices=True)
    selected_gene_ids = [model_gene_ids[index] for index in selected_indices]
    coefficients = np.asarray(
        rfecv.estimator_.named_steps["svc"].coef_, dtype=float
    ).reshape(-1)
    if len(coefficients) != len(selected_gene_ids):
        raise RuntimeError("Could not align selected genes with fitted SVC coefficients")

    # Write a complete feature table over the analysis gene universe
    model_position = {gene: index for index, gene in enumerate(model_gene_ids)}
    rank_by_gene = {
        gene: int(rfecv.ranking_[index]) for gene, index in model_position.items()
    }
    coefficient_by_gene = dict(zip(selected_gene_ids, coefficients))
    selected_gene_set = set(selected_gene_ids)
    model_gene_set = set(model_gene_ids)
    output_table = pd.DataFrame({"gene": analysis_gene_ids})
    output_table["modelled_feature"] = output_table["gene"].isin(model_gene_set).astype(int)
    output_table["selected_feature"] = output_table["gene"].isin(selected_gene_set).astype(int)
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

    # Write the actual RFECV feature counts and every fold score
    curve_table = pd.DataFrame(
        {
            "n_features": feature_counts,
            "mean_test_score": mean_scores,
            "std_test_score": std_scores,
            "selected_feature_count": (
                feature_counts == int(rfecv.n_features_)
            ).astype(int),
        }
    )
    split_score_columns = sorted(
        key
        for key in rfecv.cv_results_
        if re.fullmatch(r"split[0-9]+_test_score", key)
    )
    for column in split_score_columns:
        curve_table[column] = np.asarray(rfecv.cv_results_[column], dtype=float)
    fold_table = make_fold_table(cv, model_counts, encoded, samplesheet)

    output_prefix = arguments.output_prefix
    table_path = f"{output_prefix}_table.csv"
    curve_path = f"{output_prefix}_curve.csv"
    folds_path = f"{output_prefix}_folds.csv"
    plot_path = f"{output_prefix}_fs.png"
    model_path = f"{output_prefix}_model.pickle"
    metadata_path = f"{output_prefix}_metadata.json"
    feature_count_path = f"{output_prefix}_n_selected_features.txt"
    output_table.to_csv(table_path, index=False)
    curve_table.to_csv(curve_path, index=False)
    fold_table.to_csv(folds_path, index=False)
    Path(feature_count_path).write_text(f"{int(rfecv.n_features_)}\n", encoding="utf-8")

    # Plot the feature-selection curve against real evaluated feature counts
    os.environ.setdefault(
        "MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "zuzu-matplotlib")
    )
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots()
    axis.errorbar(
        feature_counts,
        mean_scores,
        yerr=std_scores,
        marker="o",
        capsize=3,
    )
    selected_score = mean_scores[feature_counts == int(rfecv.n_features_)][0]
    axis.scatter(
        [int(rfecv.n_features_)],
        [selected_score],
        color="red",
        zorder=3,
        label="Selected feature count",
    )
    axis.set_xlabel("Number of SVC candidate features")
    axis.set_ylabel(f"Mean cross-validation score ({arguments.scoring})")
    axis.legend()
    figure.tight_layout()
    figure.savefig(plot_path, format="png", dpi=150)
    plt.close(figure)

    resolved_parameters = {
        "kernel": "linear",
        "C": float(arguments.svc_c),
        "class_weight": resolved_class_weight,
        "scoring": arguments.scoring,
        "folds": int(arguments.folds),
        "shuffle": True,
        "seed": int(arguments.seed),
        "rfe_step": int(arguments.step),
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
        "scaling": "StandardScaler fitted inside each RFECV training fold",
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

    # Remove large fold-by-step masks added by recent scikit-learn versions
    discarded_cv_keys = []
    for key in list(rfecv.cv_results_):
        if re.fullmatch(r"split[0-9]+_(support|ranking)", key):
            discarded_cv_keys.append(key)
            del rfecv.cv_results_[key]

    # Save a model bundle that can be loaded without custom pickle classes
    model_bundle = {
        "schema_version": "zuzu_svc_model_v1",
        "rfecv": rfecv,
        "estimator": rfecv.estimator_,
        "normalization_gene_ids": normalization_gene_ids,
        "analysis_gene_ids": analysis_gene_ids,
        "model_gene_ids": model_gene_ids,
        "selected_gene_ids": selected_gene_ids,
        "class_mapping": class_mapping,
        "preprocessing": preprocessing,
        "parameters": resolved_parameters,
        "software_versions": versions,
        "discarded_cv_result_keys": discarded_cv_keys,
    }
    with open(model_path, "wb") as model_handle:
        pickle.dump(model_bundle, model_handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Save human-readable run metadata alongside the fitted model
    metadata = {
        "schema_version": "zuzu_svc_metadata_v1",
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
        },
        "selected_gene_ids": selected_gene_ids,
        "cv_feature_counts": feature_counts,
        "cv_mean_test_scores": mean_scores,
        "cv_std_test_scores": std_scores,
        "estimator_steps": list(rfecv.estimator_.named_steps),
        "pare_mode": pare_mode,
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
        },
        "software_versions": versions,
        "discarded_model_cv_result_keys": discarded_cv_keys,
    }
    with open(metadata_path, "w", encoding="utf-8") as metadata_handle:
        json.dump(json_safe(metadata), metadata_handle, indent=2, allow_nan=False)
        metadata_handle.write("\n")

    print(f"Selected {rfecv.n_features_} of {len(model_gene_ids)} modelled features")


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
    parser.add_argument("--step", required=True, type=int)
    parser.add_argument("--n-jobs", required=True, type=int)
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
