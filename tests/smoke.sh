#!/usr/bin/env bash

set -Eeuo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo_dir="$(cd "$script_dir/.." && pwd -P)"
source_sheet="$repo_dir/input/synthtest/synthsheet.csv"
source_counts="$repo_dir/input/synthtest/synthcounts.csv"
plot_dir="$script_dir/.smoke-plots"
smoke_nperms="${SMOKE_NPERMS:-1}"
smoke_synth="${SMOKE_SYNTH:-0}"
smoke_svc="${SMOKE_SVC:-0}"
smoke_max_permute_correlation="0.6"

fail() {
    printf 'Smoke test failed: %s\n' "$1" >&2
    exit 1
}

[[ -f "$repo_dir/main.nf" ]] || fail "main.nf is missing"
[[ -f "$repo_dir/nextflow.config" ]] || fail "nextflow.config is missing"
[[ -f "$source_sheet" ]] || fail "fixture sample sheet is missing"
[[ -f "$source_counts" ]] || fail "fixture count matrix is missing"
[[ "$smoke_nperms" =~ ^[1-9][0-9]*$ ]] || fail "SMOKE_NPERMS must be a positive integer"
[[ "$smoke_synth" =~ ^[01]$ ]] || fail "SMOKE_SYNTH must be either 0 or 1"
[[ "$smoke_svc" =~ ^[01]$ ]] || fail "SMOKE_SVC must be either 0 or 1"
if [[ "$smoke_svc" == 1 && "$smoke_synth" == 1 ]]; then
    fail "SMOKE_SVC=1 and SMOKE_SYNTH=1 cannot yet be combined; run the opt-in checks separately"
fi

command -v nextflow >/dev/null || fail "Nextflow is not installed"
command -v docker >/dev/null || fail "Docker is not installed"
docker info >/dev/null 2>&1 || fail "the Docker daemon is unavailable"

# KI-026 currently requires collected work paths to remain under the repository
# so they are visible inside the Docker container.
smoke_tmp_base="$repo_dir"
smoke_root="$(mktemp -d "$smoke_tmp_base/.zuzu-smoke.XXXXXX")"

compare_snapshot() {
    snapshot_source="$1"
    snapshot_copy="$2"
    snapshot_label="$3"

    if [[ ! -f "$snapshot_source" ]] || ! cmp -s "$snapshot_source" "$snapshot_copy"; then
        printf 'Smoke test failed: %s changed\n' "$snapshot_label" >&2
        return 1
    fi
}

verify_input_integrity() {
    [[ -n "${snapshot_dir:-}" && -d "$snapshot_dir" ]] || return 0

    integrity_status=0
    compare_snapshot "$source_sheet" "$snapshot_dir/source_sheet.csv" "source sample sheet" || integrity_status=1
    compare_snapshot "$source_counts" "$snapshot_dir/source_counts.csv" "source count matrix" || integrity_status=1
    compare_snapshot "$smoke_sheet" "$snapshot_dir/smoke_sheet.csv" "collision-test sample sheet" || integrity_status=1
    compare_snapshot "$smoke_counts" "$snapshot_dir/smoke_counts.tsv" "collision-test count matrix" || integrity_status=1
    return "$integrity_status"
}

cleanup() {
    smoke_status=$?
    trap - EXIT

    if ! verify_input_integrity; then
        smoke_status=1
    fi

    if [[ $smoke_status -eq 0 ]]; then
        printf 'Zuzu smoke test passed.\n'
        if [[ "${KEEP_SMOKE:-0}" == 1 ]]; then
            printf 'Temporary run retained at: %s\n' "$smoke_root"
        else
            case "$smoke_root" in
                "$smoke_tmp_base"/.zuzu-smoke.*) rm -rf -- "$smoke_root" ;;
                *) printf 'Refusing to remove unexpected path: %s\n' "$smoke_root" >&2 ;;
            esac
        fi
    else
        printf 'Smoke run retained for inspection: %s\n' "$smoke_root" >&2
    fi

    exit "$smoke_status"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

input_dir="$smoke_root/input"
launch_dir="$smoke_root/launch"
work_dir="$smoke_root/work"
output_dir="$smoke_root/output"
log_dir="$smoke_root/logs"
mkdir -p "$input_dir" "$launch_dir" "$work_dir" "$output_dir" "$log_dir"

# Deliberately reproduce the input/output basename collisions from KI-001.
# The derived fixture is also unbalanced (30 reference, 24 alternative) so the
# KI-003 regression detects accidental reuse or loss of either group size. Its
# requested reference sorts after the alternative so KI-004 cannot be masked by
# R's default lexical factor ordering.
# CLEANINPUTS must leave both source files byte-for-byte unchanged.
smoke_sheet="$input_dir/samplesheet.csv"
smoke_counts="$input_dir/countsframe_clean.csv"
smoke_reflevel="z_reference"
smoke_altlevel="a_alternative"
awk -F',' -v reference="\"$smoke_reflevel\"" -v alternative="\"$smoke_altlevel\"" '
    BEGIN { OFS = "," }
    NR > 55 { exit }
    NR == 1 { print; next }
    $3 == "\"foo\"" { $3 = reference; print; next }
    $3 == "\"nonref\"" { $3 = alternative; print; next }
    { exit 2 }
' "$source_sheet" > "$smoke_sheet" || fail "could not derive the reverse-alphabetical phenotype fixture"

# The compact tracked fixture is CSV, while CLEANINPUTS expects tab-delimited counts.
# It has no embedded commas; keep the gene column plus samples 1-54 and convert
# the result mechanically to tab-delimited input. Add four sentinel genes that
# distinguish the CPM/prevalence rule from the old raw group-mean filter.
awk -F',' 'BEGIN { OFS = "\t" }
    {
        for (i = 1; i <= 55; i++) {
            printf "%s%s", $i, (i == 55 ? ORS : OFS)
        }
    }
    END {
        printf "smoke_filter_below_cpm"
        for (i = 2; i <= 55; i++) printf "%s6", OFS
        printf "%s", ORS

        printf "smoke_filter_prevalent"
        for (i = 2; i <= 55; i++) printf "%s20", OFS
        printf "%s", ORS

        printf "smoke_filter_group_specific"
        for (i = 2; i <= 55; i++) printf "%s%d", OFS, (i >= 32 ? 20 : 0)
        printf "%s", ORS

        printf "smoke_filter_sparse"
        for (i = 2; i <= 55; i++) printf "%s%d", OFS, (i <= 24 ? 1000 : 0)
        printf "%s", ORS

        # These rows are discarded before analysis and make the high-filtering
        # warning deterministic even if every tracked fixture gene is retained.
        for (sentinel = 1; sentinel <= 1001; sentinel++) {
            printf "smoke_filter_bulk_low_%04d", sentinel
            for (i = 2; i <= 55; i++) printf "%s0", OFS
            printf "%s", ORS
        }
    }
' "$source_counts" > "$smoke_counts"

snapshot_dir="$smoke_root/snapshots"
mkdir -p "$snapshot_dir"
cp "$source_sheet" "$snapshot_dir/source_sheet.csv"
cp "$source_counts" "$snapshot_dir/source_counts.csv"
cp "$smoke_sheet" "$snapshot_dir/smoke_sheet.csv"
cp "$smoke_counts" "$snapshot_dir/smoke_counts.tsv"

synth_args=(--synthstep false)
if [[ "$smoke_synth" == 1 ]]; then
    synth_args=(--synthstep true --synthdeg_pct 0.1)
fi

svc_args=(--mlstep false)
if [[ "$smoke_svc" == 1 ]]; then
    svc_args=(--mlstep true)
fi

# The relaxed test-only anchor threshold gives the main smoke run stable,
# nonempty anchor sets. Focused checks below cover the zero- and low-anchor cases.
(
    cd "$launch_dir"
    NXF_ANSI_LOG=false nextflow \
        -log "$log_dir/nextflow.debug.log" \
        -C "$repo_dir/nextflow.config" \
        run "$repo_dir/main.nf" \
        -profile local \
        --samplesheet "$smoke_sheet" \
        --countsframe "$smoke_counts" \
        --reflevel "$smoke_reflevel" \
        --outdir "$output_dir" \
        --nperms "$smoke_nperms" \
        --truecutoff 0.05 \
        "${synth_args[@]}" \
        "${svc_args[@]}" \
        --pareml true \
        -w "$work_dir" \
        -with-trace "$log_dir/trace.tsv"
) 2>&1 | tee "$log_dir/console.log"

verify_input_integrity || fail "input integrity check failed"

require_file() {
    [[ -s "$1" ]] || fail "missing or empty artifact: $1"
}

require_pdf() {
    require_file "$1"
    LC_ALL=C head -c 5 "$1" | grep -q '^%PDF-' || fail "not a valid PDF: $1"
}

require_png() {
    require_file "$1"
    png_signature="$(od -An -tx1 -N8 "$1" | tr -d '[:space:]')"
    [[ "$png_signature" == "89504e470d0a1a0a" ]] || fail "not a valid PNG: $1"
}

csv_value() {
    csv_file="$1"
    wanted_row="$2"
    wanted_column="$3"

    awk -F',' -v wanted_row="$wanted_row" -v wanted_column="$wanted_column" '
        function unquote(value) {
            sub(/^"/, "", value)
            sub(/"$/, "", value)
            return value
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if (unquote($i) == wanted_column) {
                    value_column = i
                }
            }
            if (!value_column) {
                exit 2
            }
            next
        }
        unquote($1) == wanted_row {
            print unquote($(value_column))
            found_row = 1
            exit
        }
        END {
            if (!found_row) {
                exit 3
            }
        }
    ' "$csv_file"
}

locate_clean_outputs() {
    clean_counts=""
    clean_counts_matches=0
    while IFS= read -r candidate_counts; do
        clean_counts="$candidate_counts"
        ((clean_counts_matches += 1))
    done < <(find "$work_dir" -mindepth 3 -maxdepth 3 -type f -name 'countsframe_clean.csv' -print)

    [[ "$clean_counts_matches" -eq 1 ]] ||
        fail "expected one cleaned count matrix, found $clean_counts_matches"
    clean_task_dir="${clean_counts%/*}"
    clean_sheet="$clean_task_dir/samplesheet.csv"
    require_file "$clean_counts"
    require_file "$clean_sheet"
}

verify_common_filter() {
    prevalent_sample1="$(csv_value "$clean_counts" "smoke_filter_prevalent" "sample1")" ||
        fail "the prevalent filtering sentinel was removed"
    prevalent_sample54="$(csv_value "$clean_counts" "smoke_filter_prevalent" "sample54")" ||
        fail "the prevalent filtering sentinel was removed"
    [[ "$prevalent_sample1:$prevalent_sample54" == "20:20" ]] ||
        fail "filtering changed counts for the prevalent sentinel"

    group_specific_sample1="$(csv_value "$clean_counts" "smoke_filter_group_specific" "sample1")" ||
        fail "the group-specific filtering sentinel was removed"
    group_specific_sample31="$(csv_value "$clean_counts" "smoke_filter_group_specific" "sample31")" ||
        fail "the group-specific filtering sentinel was removed"
    [[ "$group_specific_sample1:$group_specific_sample31" == "0:20" ]] ||
        fail "filtering changed counts for the group-specific sentinel"

    if csv_value "$clean_counts" "smoke_filter_below_cpm" "sample1" >/dev/null 2>&1; then
        fail "the below-CPM filtering sentinel was retained"
    fi
    if csv_value "$clean_counts" "smoke_filter_sparse" "sample1" >/dev/null 2>&1; then
        fail "the low-prevalence filtering sentinel was retained"
    fi
}

verify_filter_warning() {
    warning_prefix="WARNING [CLEANINPUTS]: the common CPM/prevalence filter removed"
    grep -Fq "$warning_prefix" "$log_dir/console.log" ||
        fail "the high-filtering warning was not reported in the console"
    grep -Fq "exceeding the 50% warning threshold" "$log_dir/console.log" ||
        fail "the high-filtering warning did not report its threshold"
}

verify_reference_direction() {
    deseq_basic_table=""
    deseq_basic_matches=0
    while IFS= read -r dds_file; do
        candidate_dir="${dds_file%/*}"
        candidate_table="$candidate_dir/deseq_table.csv"
        [[ -f "$candidate_table" ]] || continue
        deseq_basic_table="$candidate_table"
        ((deseq_basic_matches += 1))
    done < <(find "$work_dir" -mindepth 3 -maxdepth 3 -type f -name 'dds_deg.RData' -print)

    [[ "$deseq_basic_matches" -eq 1 ]] ||
        fail "expected one DESeq2 basic result table, found $deseq_basic_matches"

    edger_basic_table=""
    edger_basic_matches=0
    while IFS= read -r candidate_table; do
        candidate_dir="${candidate_table%/*}"
        [[ ! -e "$candidate_dir/nDEGs.RData" ]] || continue
        edger_basic_table="$candidate_table"
        ((edger_basic_matches += 1))
    done < <(find "$work_dir" -mindepth 3 -maxdepth 3 -type f -name 'edger_table.csv' -print)

    [[ "$edger_basic_matches" -eq 1 ]] ||
        fail "expected one edgeR basic result table, found $edger_basic_matches"

    # g2 is substantially higher in the alternative samples in this fixture.
    # Both methods must therefore report a positive alternative-vs-reference effect.
    deseq_effect="$(csv_value "$deseq_basic_table" "g2" "log2FoldChange")" ||
        fail "could not read the DESeq2 g2 effect"
    edger_effect="$(csv_value "$edger_basic_table" "g2" "logFC")" ||
        fail "could not read the edgeR g2 effect"

    awk -v effect="$deseq_effect" 'BEGIN { exit !(effect > 0) }' ||
        fail "DESeq2 did not report g2 as higher in the alternative group: $deseq_effect"
    awk -v effect="$edger_effect" 'BEGIN { exit !(effect > 0) }' ||
        fail "edgeR did not report g2 as higher in the alternative group: $edger_effect"
}

permutation_sheet_stats() {
    original_sheet="$1"
    permuted_sheet="$2"

    awk -F',' -v expected_ref="$smoke_reflevel" -v expected_alt="$smoke_altlevel" '
        function unquote(value) {
            sub(/^"/, "", value)
            sub(/"$/, "", value)
            return value
        }
        function abort_with(code) {
            failure = code
            exit code
        }
        FNR == 1 {
            file_number += 1
            for (i = 1; i <= NF; i++) {
                heading = unquote($i)
                if (heading == "sample") {
                    sample_column[file_number] = i
                } else if (heading == "phenotype") {
                    phenotype_column[file_number] = i
                }
            }
            if (!sample_column[file_number] || !phenotype_column[file_number]) {
                abort_with(2)
            }
            next
        }
        file_number == 1 {
            sample = unquote($(sample_column[file_number]))
            phenotype = unquote($(phenotype_column[file_number]))
            if (sample in original_phenotype) {
                abort_with(3)
            }
            if (phenotype != expected_ref && phenotype != expected_alt) {
                abort_with(4)
            }
            original_rows += 1
            original_order[original_rows] = sample
            original_phenotype[sample] = phenotype
            original_ref += (phenotype == expected_ref)
            original_alt += (phenotype == expected_alt)
            next
        }
        file_number == 2 {
            sample = unquote($(sample_column[file_number]))
            phenotype = unquote($(phenotype_column[file_number]))
            if (!(sample in original_phenotype) || sample in seen_permuted) {
                abort_with(5)
            }
            permuted_rows += 1
            if (sample != original_order[permuted_rows]) {
                abort_with(6)
            }
            if (phenotype != expected_ref && phenotype != expected_alt) {
                abort_with(7)
            }
            seen_permuted[sample] = 1
            permuted_ref += (phenotype == expected_ref)
            permuted_alt += (phenotype == expected_alt)
            both_alt += (original_phenotype[sample] == expected_alt && phenotype == expected_alt)
            next
        }
        END {
            if (failure) {
                exit failure
            }
            if (file_number != 2 || original_rows != permuted_rows) {
                exit 8
            }
            if (original_ref != permuted_ref || original_alt != permuted_alt) {
                exit 9
            }
            denominator = sqrt(original_alt * original_ref * permuted_alt * permuted_ref)
            if (denominator == 0) {
                exit 10
            }
            correlation = (original_rows * both_alt - original_alt * permuted_alt) / denominator
            printf "%.12g %d %d\n", correlation, permuted_ref, permuted_alt
        }
    ' "$original_sheet" "$permuted_sheet"
}

verify_null_permutations() {
    unexpected_permuted_counts="$(find "$work_dir" -name 'countsframe_perm.csv' -print -quit)"
    [[ -z "$unexpected_permuted_counts" ]] ||
        fail "null permutation created a replacement count matrix: $unexpected_permuted_counts"

    permutation_count=0
    permutation_replicate_counts=()
    while IFS= read -r metadata_file; do
        task_dir="${metadata_file%/*}"
        permuted_sheet="$task_dir/samplesheet_perm.csv"
        forwarded_counts="$task_dir/countsframe_clean.csv"
        require_file "$permuted_sheet"
        require_file "$forwarded_counts"
        cmp -s "$clean_counts" "$forwarded_counts" ||
            fail "null permutation changed the cleaned count matrix: $forwarded_counts"

        sheet_stats="$(permutation_sheet_stats "$clean_sheet" "$permuted_sheet")" ||
            fail "invalid null-permutation sample sheet: $permuted_sheet"
        read -r observed_correlation observed_ref_count observed_alt_count <<< "$sheet_stats"
        [[ "$observed_ref_count:$observed_alt_count" == "30:24" ]] ||
            fail "null permutation changed phenotype counts: $observed_ref_count ref, $observed_alt_count alt"
        awk -v value="$observed_correlation" -v maximum="$smoke_max_permute_correlation" '
            BEGIN { exit !(value < maximum && value > -maximum) }
        ' || fail "null permutation correlation exceeds the constraint: $observed_correlation"

        metadata_values="$(awk -F',' '
            NR == 1 {
                for (i = 1; i <= NF; i++) {
                    gsub(/^"|"$/, "", $i)
                    column[$i] = i
                }
                required[1] = "replicate"
                required[2] = "correlation"
                required[3] = "max_abs_correlation"
                required[4] = "reference"
                required[5] = "alternative"
                required[6] = "reference_count"
                required[7] = "alternative_count"
                required[8] = "alternative_overlap"
                for (i = 1; i <= 8; i++) {
                    if (!(required[i] in column)) {
                        exit 2
                    }
                }
                next
            }
            NR == 2 {
                for (i = 1; i <= NF; i++) {
                    gsub(/^"|"$/, "", $i)
                }
                printf "%s %s %s %s %s %s %s %s\n",
                       $(column["replicate"]),
                       $(column["correlation"]),
                       $(column["max_abs_correlation"]),
                       $(column["reference"]),
                       $(column["alternative"]),
                       $(column["reference_count"]),
                       $(column["alternative_count"]),
                       $(column["alternative_overlap"])
                next
            }
            END {
                if (NR != 2) {
                    exit 3
                }
            }
        ' "$metadata_file")" || fail "invalid permutation metadata: $metadata_file"

        read -r replicate_id metadata_correlation metadata_maximum metadata_ref metadata_alt \
            metadata_ref_count metadata_alt_count metadata_overlap <<< "$metadata_values"
        [[ "$replicate_id" =~ ^[1-9][0-9]*$ ]] || fail "invalid permutation replicate: $replicate_id"
        [[ "$metadata_maximum" == "$smoke_max_permute_correlation" ]] ||
            fail "unexpected permutation correlation limit: $metadata_maximum"
        [[ "$metadata_ref:$metadata_alt" == "$smoke_reflevel:$smoke_altlevel" ]] ||
            fail "unexpected permutation phenotype labels: $metadata_ref, $metadata_alt"
        [[ "$metadata_ref_count:$metadata_alt_count" == "30:24" ]] ||
            fail "invalid permutation metadata counts: $metadata_ref_count ref, $metadata_alt_count alt"
        [[ "$metadata_overlap" =~ ^[0-9]+$ ]] || fail "invalid alternative overlap: $metadata_overlap"
        awk -v observed="$observed_correlation" -v recorded="$metadata_correlation" '
            BEGIN {
                difference = observed - recorded
                if (difference < 0) difference = -difference
                exit !(difference < 1e-10)
            }
        ' || fail "permutation metadata correlation does not match its sample sheet"

        permutation_replicate_counts[$replicate_id]=$((
            ${permutation_replicate_counts[$replicate_id]:-0} + 1
        ))
        ((permutation_count += 1))
    done < <(find "$work_dir" -mindepth 3 -maxdepth 3 -type f -name 'permutation_metadata.csv' -print | sort)

    [[ "$permutation_count" -eq "$smoke_nperms" ]] ||
        fail "expected $smoke_nperms null permutations, found $permutation_count"
    for ((replicate_id = 1; replicate_id <= smoke_nperms; replicate_id++)); do
        [[ "${permutation_replicate_counts[$replicate_id]:-0}" -eq 1 ]] ||
            fail "expected one null permutation for replicate $replicate_id"
    done
}

verify_permutation_trace_keys() {
    for process_name in DATA_PERMUTE DESEQ_PERMUTE EDGER_PERMUTE WILCOXON_PERMUTE; do
        for ((replicate_id = 1; replicate_id <= smoke_nperms; replicate_id++)); do
            expected_tag="$process_name (replicate=$replicate_id)"
            trace_matches="$(awk -F '\t' -v expected="$expected_tag" '
                NR == 1 {
                    for (i = 1; i <= NF; i++) {
                        if ($i == "name") {
                            name_column = i
                        } else if ($i == "status") {
                            status_column = i
                        }
                    }
                    if (!name_column || !status_column) {
                        exit 2
                    }
                    next
                }
                index($(name_column), expected) && $(status_column) == "COMPLETED" {
                    matches += 1
                }
                END {
                    print matches + 0
                }
            ' "$log_dir/trace.tsv")" || fail "could not parse permutation task keys from the trace"

            [[ "$trace_matches" -eq 1 ]] ||
                fail "expected one completed task for $expected_tag, found $trace_matches"
        done
    done
}

validate_quasi_metric_rows() {
    metrics_file="$1"
    expected_rows="$2"
    require_method="$3"

    awk -F',' -v expected_rows="$expected_rows" -v require_method="$require_method" '
        function absolute(value) {
            return value < 0 ? -value : value
        }
        function is_integer(value) {
            return value ~ /^[0-9]+$/
        }
        function is_number(value) {
            return value ~ /^-?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                gsub(/^"|"$/, "", $i)
                column[$i] = i
            }
            required["power"] = 1
            required["FDP"] = 1
            required["youden_index"] = 1
            required["true_positives"] = 1
            required["false_positives"] = 1
            required["true_negatives"] = 1
            required["false_negatives"] = 1
            required["n_tested"] = 1
            required["n_called"] = 1
            required["n_truth"] = 1
            required["n_background"] = 1
            required["metric_status"] = 1
            required["replicate"] = 1
            required["refnum"] = 1
            required["altnum"] = 1
            required["totalnum"] = 1
            for (name in required) {
                if (!(name in column)) exit 2
            }
            if (("normAUC" in column) || (require_method && !("method" in column))) exit 3
            next
        }
        {
            data_rows += 1
            for (i = 1; i <= NF; i++) gsub(/^"|"$/, "", $i)

            power = $(column["power"])
            FDP = $(column["FDP"])
            youden = $(column["youden_index"])
            TP = $(column["true_positives"])
            FP = $(column["false_positives"])
            TN = $(column["true_negatives"])
            FN = $(column["false_negatives"])
            n_tested = $(column["n_tested"])
            n_called = $(column["n_called"])
            n_truth = $(column["n_truth"])
            n_background = $(column["n_background"])
            status = $(column["metric_status"])
            replicate = $(column["replicate"])
            refnum = $(column["refnum"])
            altnum = $(column["altnum"])
            totalnum = $(column["totalnum"])

            if (!is_number(FDP)) exit 4
            if (!is_integer(TP) || !is_integer(FP) || !is_integer(TN) || !is_integer(FN) ||
                !is_integer(n_tested) || !is_integer(n_called) || !is_integer(n_truth) ||
                !is_integer(n_background) || !is_integer(replicate) || !is_integer(refnum) ||
                !is_integer(altnum) || !is_integer(totalnum)) exit 5
            if (n_tested == 0 || n_tested != TP + FP + TN + FN ||
                n_called != TP + FP || n_truth != TP + FN ||
                n_background != FP + TN || n_tested != n_truth + n_background) exit 6
            if (totalnum != refnum + altnum) exit 7

            expected_FDP = n_called == 0 ? 0 : FP / n_called
            if (absolute(FDP - expected_FDP) > 1e-9 || FDP < 0 || FDP > 1) exit 8

            if (n_truth == 0) {
                if (power != "NA" || youden != "NA") exit 9
            } else {
                if (!is_number(power)) exit 10
                expected_power = TP / n_truth
                if (absolute(power - expected_power) > 1e-9 || power < 0 || power > 1) exit 11

                if (n_background == 0) {
                    if (youden != "NA") exit 12
                } else {
                    if (!is_number(youden)) exit 13
                    expected_youden = expected_power - (FP / n_background)
                    if (absolute(youden - expected_youden) > 1e-9 ||
                        youden < -1 || youden > 1) exit 14
                }
            }

            if (n_truth == 0) {
                expected_status = "no_anchors"
            } else if (n_background == 0) {
                expected_status = "youden_undefined_no_background"
            } else if (n_truth <= 10) {
                expected_status = "low_anchor_count"
            } else {
                expected_status = n_called == 0 ? "no_calls" : "ok"
            }
            if (status != expected_status) exit 15

            method = require_method ? $(column["method"]) : ""
            if (require_method) {
                if (method != "DESeq2" && method != "edgeR" && method != "Wilcoxon") exit 16
                key = method SUBSEP replicate SUBSEP refnum SUBSEP altnum
                if (seen[key]++) exit 17
                method_count[method] += 1
            }
            printf "%d %d %s\n", refnum, altnum, method
        }
        END {
            if (data_rows != expected_rows) exit 18
            if (require_method) {
                expected_per_method = expected_rows / 3
                if (method_count["DESeq2"] != expected_per_method ||
                    method_count["edgeR"] != expected_per_method ||
                    method_count["Wilcoxon"] != expected_per_method) exit 19
            }
        }
    ' "$metrics_file"
}

# Verify the zero-, low-, and adequate-anchor metric conventions directly
verify_quasi_anchor_contract() {
    contract_file="$smoke_root/quasi_anchor_contract.csv"
    invalid_contract_file="$smoke_root/quasi_anchor_contract_invalid.csv"

    printf '%s\n' \
        'power,FDP,youden_index,true_positives,false_positives,true_negatives,false_negatives,n_tested,n_called,n_truth,n_background,metric_status,replicate,refnum,altnum,totalnum' \
        'NA,0,NA,0,0,20,0,20,0,0,20,no_anchors,1,5,5,10' \
        'NA,1,NA,0,2,18,0,20,2,0,20,no_anchors,2,5,5,10' \
        '1,0.5,0.947368421052632,1,1,18,0,20,2,1,19,low_anchor_count,3,5,5,10' \
        '0,0,0,0,0,10,10,20,0,10,10,low_anchor_count,4,5,5,10' \
        '0,0,0,0,0,9,11,20,0,11,9,no_calls,5,5,5,10' \
        '0.4,0,NA,2,0,0,3,5,2,5,0,youden_undefined_no_background,6,5,5,10' \
        > "$contract_file"

    validate_quasi_metric_rows "$contract_file" 6 0 >/dev/null ||
        fail "quasi anchor-count metric contract is invalid"

    printf '%s\n' \
        'power,FDP,youden_index,true_positives,false_positives,true_negatives,false_negatives,n_tested,n_called,n_truth,n_background,metric_status,replicate,refnum,altnum,totalnum' \
        '0,0,0,0,0,20,0,20,0,0,20,no_anchors,1,5,5,10' \
        > "$invalid_contract_file"

    if validate_quasi_metric_rows "$invalid_contract_file" 1 0 >/dev/null 2>&1; then
        fail "quasi validator accepted defined recovery metrics with no anchors"
    fi

    for quasi_module in deseq_quasi.nf edger_quasi.nf wilcoxon_quasi.nf; do
        module_path="$repo_dir/modules/$quasi_module"
        grep -q '"no_anchors"' "$module_path" ||
            fail "$quasi_module does not report the no-anchor state"
        grep -q '"low_anchor_count"' "$module_path" ||
            fail "$quasi_module does not report the low-anchor state"
        if grep -Eq 'length\(truedegs\)[[:space:]]*>?=[[:space:]]*10|length\(truedegs\)[[:space:]]*>[[:space:]]*10|insufficient_truth' "$module_path"; then
            fail "$quasi_module still censors low-anchor datasets"
        fi
    done
}

verify_quasi_metrics() {
    expected_result_count=$((smoke_nperms * 3 * 3))
    expected_pair_count=$((smoke_nperms * 3))
    quasi_result_count=0
    small_pair_count=0
    middle_pair_count=0
    full_pair_count=0

    while IFS= read -r quasi_result; do
        if ! grep -q '"refnum".*"altnum".*"totalnum"' "$quasi_result"; then
            continue
        fi
        if head -n 1 "$quasi_result" | grep -q '"depth"'; then
            continue
        fi

        metric_values="$(validate_quasi_metric_rows "$quasi_result" 1 0)" ||
            fail "invalid quasi metric row: $quasi_result"
        read -r observed_ref observed_alt _ <<< "$metric_values"
        ((quasi_result_count += 1))
        case "$observed_ref:$observed_alt" in
            11:9) ((small_pair_count += 1)) ;;
            20:16) ((middle_pair_count += 1)) ;;
            30:24) ((full_pair_count += 1)) ;;
            *) fail "unexpected quasi sample-size pair: $observed_ref ref, $observed_alt alt" ;;
        esac
    done < <(find "$work_dir" -mindepth 3 -maxdepth 3 -type f -name 'outframe.csv' -print | sort)

    [[ "$quasi_result_count" -eq "$expected_result_count" ]] ||
        fail "expected $expected_result_count quasi result rows, found $quasi_result_count"
    [[ "$small_pair_count" -eq "$expected_pair_count" ]] ||
        fail "expected $expected_pair_count quasi results for 11 ref + 9 alt, found $small_pair_count"
    [[ "$middle_pair_count" -eq "$expected_pair_count" ]] ||
        fail "expected $expected_pair_count quasi results for 20 ref + 16 alt, found $middle_pair_count"
    [[ "$full_pair_count" -eq "$expected_pair_count" ]] ||
        fail "expected $expected_pair_count quasi results for 30 ref + 24 alt, found $full_pair_count"

    validate_quasi_metric_rows "$output_dir/quasi_results.csv" "$expected_result_count" 1 >/dev/null ||
        fail "invalid combined quasi metrics: $output_dir/quasi_results.csv"
}

verify_fullsynth_metrics() {
    expected_result_count=$((smoke_nperms * 3 * 3 * 3))
    expected_summary_count=$((3 * 3 * 3 * 3))

    awk -F',' -v expected_rows="$expected_result_count" -v expected_replicates="$smoke_nperms" '
        function absolute(value) { return value < 0 ? -value : value }
        function is_integer(value) { return value ~ /^[0-9]+$/ }
        function is_number(value) {
            return value ~ /^-?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                gsub(/^"|"$/, "", $i)
                column[$i] = i
            }
            split("method power FDP youden_index true_positives false_positives true_negatives false_negatives n_tested n_called n_truth n_background metric_status replicate refnum altnum totalnum depth", names, " ")
            for (i in names) if (!(names[i] in column)) exit 2
            next
        }
        {
            rows += 1
            for (i = 1; i <= NF; i++) gsub(/^"|"$/, "", $i)
            method = $(column["method"])
            power = $(column["power"])
            FDP = $(column["FDP"])
            youden = $(column["youden_index"])
            TP = $(column["true_positives"])
            FP = $(column["false_positives"])
            TN = $(column["true_negatives"])
            FN = $(column["false_negatives"])
            n_tested = $(column["n_tested"])
            n_called = $(column["n_called"])
            n_truth = $(column["n_truth"])
            n_background = $(column["n_background"])
            status = $(column["metric_status"])
            replicate = $(column["replicate"])
            refnum = $(column["refnum"])
            altnum = $(column["altnum"])
            totalnum = $(column["totalnum"])
            depth = $(column["depth"])

            if (method != "DESeq2" && method != "edgeR" && method != "Wilcoxon") exit 3
            if (!is_number(power) || !is_number(FDP) || !is_number(youden) ||
                !is_number(depth)) exit 4
            if (!is_integer(TP) || !is_integer(FP) || !is_integer(TN) || !is_integer(FN) ||
                !is_integer(n_tested) || !is_integer(n_called) || !is_integer(n_truth) ||
                !is_integer(n_background) || !is_integer(replicate) ||
                !is_integer(refnum) || !is_integer(altnum) || !is_integer(totalnum)) exit 5
            if (n_truth == 0 || n_background == 0 ||
                n_tested != TP + FP + TN + FN || n_called != TP + FP ||
                n_truth != TP + FN || n_background != FP + TN ||
                n_tested != n_truth + n_background || totalnum != refnum + altnum ||
                replicate < 1 || replicate > expected_replicates) exit 6

            pair = refnum ":" altnum
            if (pair != "11:9" && pair != "20:16" && pair != "30:24") exit 7
            if (depth != 100000 && depth != 1000000 && depth != 10000000) exit 8

            expected_status = n_called == 0 ? "no_calls" : "ok"
            expected_power = TP / n_truth
            expected_FDP = n_called == 0 ? 0 : FP / n_called
            expected_youden = expected_power - (FP / n_background)
            if (status != expected_status || absolute(power - expected_power) > 1e-9 ||
                absolute(FDP - expected_FDP) > 1e-9 ||
                absolute(youden - expected_youden) > 1e-9) exit 9

            key = method SUBSEP replicate SUBSEP refnum SUBSEP altnum SUBSEP depth
            if (seen[key]++) exit 10
            dataset_key = replicate SUBSEP refnum SUBSEP altnum SUBSEP depth
            if (dataset_key in dataset_tested &&
                (dataset_tested[dataset_key] != n_tested ||
                 dataset_truth[dataset_key] != n_truth ||
                 dataset_background[dataset_key] != n_background)) exit 11
            dataset_tested[dataset_key] = n_tested
            dataset_truth[dataset_key] = n_truth
            dataset_background[dataset_key] = n_background
            method_count[method] += 1
            pair_count[pair] += 1
            method_pair_count[method SUBSEP pair] += 1
        }
        END {
            if (rows != expected_rows) exit 12
            expected_per_method = expected_rows / 3
            if (method_count["DESeq2"] != expected_per_method ||
                method_count["edgeR"] != expected_per_method ||
                method_count["Wilcoxon"] != expected_per_method) exit 13

            expected_per_pair = expected_rows / 3
            expected_per_method_pair = expected_rows / 9
            split("11:9 20:16 30:24", expected_pairs, " ")
            split("DESeq2 edgeR Wilcoxon", expected_methods, " ")
            for (i in expected_pairs) {
                pair = expected_pairs[i]
                if (pair_count[pair] != expected_per_pair) exit 14
                for (j in expected_methods) {
                    method = expected_methods[j]
                    if (method_pair_count[method SUBSEP pair] != expected_per_method_pair) exit 15
                }
            }
        }
    ' "$output_dir/synth_results.csv" || fail "invalid synthetic per-replicate metrics"

    awk -F',' -v expected_rows="$expected_summary_count" -v expected_replicates="$smoke_nperms" '
        function is_integer(value) { return value ~ /^[0-9]+$/ }
        function is_number(value) {
            return value ~ /^-?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                gsub(/^"|"$/, "", $i)
                column[$i] = i
            }
            split("method refnum altnum totalnum depth variable mean n_valid n_total n_no_calls", names, " ")
            for (i in names) if (!(names[i] in column)) exit 2
            next
        }
        {
            rows += 1
            for (i = 1; i <= NF; i++) gsub(/^"|"$/, "", $i)
            method = $(column["method"])
            refnum = $(column["refnum"])
            altnum = $(column["altnum"])
            totalnum = $(column["totalnum"])
            depth = $(column["depth"])
            variable = $(column["variable"])
            mean = $(column["mean"])
            n_valid = $(column["n_valid"])
            n_total = $(column["n_total"])
            n_no_calls = $(column["n_no_calls"])

            if (method != "DESeq2" && method != "edgeR" && method != "Wilcoxon") exit 3
            if (variable != "power" && variable != "FDP" && variable != "youden_index") exit 4
            if (!is_number(mean) || !is_integer(n_valid) || !is_integer(n_total) ||
                !is_integer(n_no_calls) || !is_integer(refnum) || !is_integer(altnum) ||
                !is_integer(totalnum) || !is_number(depth)) exit 5
            if (n_valid != expected_replicates || n_total != expected_replicates ||
                n_no_calls > n_total || totalnum != refnum + altnum) exit 6
            if ((variable == "power" || variable == "FDP") && (mean < 0 || mean > 1)) exit 7
            if (variable == "youden_index" && (mean < -1 || mean > 1)) exit 8

            pair = refnum ":" altnum
            if (pair != "11:9" && pair != "20:16" && pair != "30:24") exit 9
            if (depth != 100000 && depth != 1000000 && depth != 10000000) exit 10

            key = method SUBSEP refnum SUBSEP altnum SUBSEP depth SUBSEP variable
            if (seen[key]++) exit 11
            group_key = method SUBSEP refnum SUBSEP altnum SUBSEP depth
            if (group_key in group_no_calls && group_no_calls[group_key] != n_no_calls) exit 12
            group_no_calls[group_key] = n_no_calls
            pair_count[pair] += 1
            method_pair_count[method SUBSEP pair] += 1
        }
        END {
            if (rows != expected_rows) exit 13

            expected_per_pair = expected_rows / 3
            expected_per_method_pair = expected_rows / 9
            split("11:9 20:16 30:24", expected_pairs, " ")
            split("DESeq2 edgeR Wilcoxon", expected_methods, " ")
            for (i in expected_pairs) {
                pair = expected_pairs[i]
                if (pair_count[pair] != expected_per_pair) exit 14
                for (j in expected_methods) {
                    method = expected_methods[j]
                    if (method_pair_count[method SUBSEP pair] != expected_per_method_pair) exit 15
                }
            }
        }
    ' "$output_dir/synth_summary.csv" || fail "invalid synthetic metric summary"
}

# Verify that each synthetic dataset retains the exact requested group sizes
# and that its count columns remain aligned with its sample sheet.
verify_fullsynth_datasets() {
    expected_dataset_count=$((smoke_nperms * 3 * 3))
    expected_pair_count=$((smoke_nperms * 3))
    dataset_count=0
    small_pair_count=0
    middle_pair_count=0
    full_pair_count=0

    while IFS= read -r synth_sheet; do
        synth_dir="$(dirname "$synth_sheet")"
        synth_counts="$synth_dir/synthcounts.csv"
        [[ -f "$synth_counts" ]] || fail "synthetic sample sheet has no matching count matrix: $synth_sheet"

        dataset_sizes="$(awk -F',' -v reflevel="$smoke_reflevel" '
            function clean(value) {
                sub(/\r$/, "", value)
                sub(/^"/, "", value)
                sub(/"$/, "", value)
                return value
            }
            FILENAME == ARGV[1] {
                if (FNR == 1) {
                    if (NF != 2 || clean($1) != "sample" || clean($2) != "phenotype") exit 2
                    sheet_header = 1
                    next
                }

                sample = clean($1)
                phenotype = clean($2)
                if (NF != 2 || sample == "" || phenotype == "") exit 3
                if (sample_seen[sample]++) exit 4
                sheet_sample[++sample_count] = sample
                if (phenotype == reflevel) {
                    ref_count += 1
                } else {
                    if (altlevel == "") altlevel = phenotype
                    if (phenotype != altlevel) exit 5
                    alt_count += 1
                }
                next
            }
            FILENAME == ARGV[2] {
                if (FNR == 1) {
                    if (NF != sample_count + 1 || clean($1) != "") exit 6
                    for (i = 2; i <= NF; i++) {
                        if (clean($i) != sheet_sample[i - 1]) exit 7
                    }
                    counts_header = 1
                    next
                }

                gene = clean($1)
                if (NF != sample_count + 1 || gene == "" || gene_seen[gene]++) exit 8
                gene_count += 1
            }
            END {
                if (!sheet_header || !counts_header || sample_count == 0 || gene_count == 0 ||
                    ref_count == 0 || alt_count == 0 || ref_count + alt_count != sample_count) exit 9
                print ref_count, alt_count, sample_count
            }
        ' "$synth_sheet" "$synth_counts")" ||
            fail "invalid synthetic sample/count pair: $synth_dir"

        read -r observed_ref observed_alt observed_total <<< "$dataset_sizes"
        [[ "$observed_total" -eq $((observed_ref + observed_alt)) ]] ||
            fail "inconsistent synthetic total in $synth_dir"
        ((dataset_count += 1))
        case "$observed_ref:$observed_alt" in
            11:9) ((small_pair_count += 1)) ;;
            20:16) ((middle_pair_count += 1)) ;;
            30:24) ((full_pair_count += 1)) ;;
            *) fail "unexpected synthetic sample-size pair: $observed_ref ref, $observed_alt alt" ;;
        esac
    done < <(find "$work_dir" -mindepth 3 -maxdepth 3 -type f -name 'synthsheet.csv' -print | sort)

    synth_count_count="$(find "$work_dir" -mindepth 3 -maxdepth 3 -type f -name 'synthcounts.csv' -print | awk 'END { print NR + 0 }')"
    [[ "$dataset_count" -eq "$expected_dataset_count" ]] ||
        fail "expected $expected_dataset_count generated synthetic datasets, found $dataset_count"
    [[ "$synth_count_count" -eq "$expected_dataset_count" ]] ||
        fail "expected $expected_dataset_count generated synthetic count matrices, found $synth_count_count"
    [[ "$small_pair_count" -eq "$expected_pair_count" ]] ||
        fail "expected $expected_pair_count synthetic datasets for 11 ref + 9 alt, found $small_pair_count"
    [[ "$middle_pair_count" -eq "$expected_pair_count" ]] ||
        fail "expected $expected_pair_count synthetic datasets for 20 ref + 16 alt, found $middle_pair_count"
    [[ "$full_pair_count" -eq "$expected_pair_count" ]] ||
        fail "expected $expected_pair_count synthetic datasets for 30 ref + 24 alt, found $full_pair_count"

    completed_tasks="$(awk -F '\t' '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "name") name_column = i
                if ($i == "status") status_column = i
            }
            if (!name_column || !status_column) exit 2
            next
        }
        index($(name_column), "DATA_FULLSYNTH (") && $(status_column) == "COMPLETED" {
            matches += 1
        }
        END { print matches + 0 }
    ' "$log_dir/trace.tsv")" || fail "could not parse synthetic generation tasks from the trace"
    [[ "$completed_tasks" -eq "$expected_dataset_count" ]] ||
        fail "expected $expected_dataset_count completed DATA_FULLSYNTH tasks, found $completed_tasks"
}

verify_quasi_trace_keys() {
    for process_name in DESEQ_QUASI EDGER_QUASI WILCOXON_QUASI; do
        for ((replicate_id = 1; replicate_id <= smoke_nperms; replicate_id++)); do
            for size_pair in "11 9" "20 16" "30 24"; do
                read -r expected_ref expected_alt <<< "$size_pair"
                expected_tag="$process_name (replicate=$replicate_id,ref=$expected_ref,alt=$expected_alt)"
                trace_matches="$(awk -F '\t' -v expected="$expected_tag" '
                    NR == 1 {
                        for (i = 1; i <= NF; i++) {
                            if ($i == "name") {
                                name_column = i
                            } else if ($i == "status") {
                                status_column = i
                            }
                        }
                        if (!name_column || !status_column) {
                            exit 2
                        }
                        next
                    }
                    index($(name_column), expected) && $(status_column) == "COMPLETED" {
                        matches += 1
                    }
                    END {
                        print matches + 0
                    }
                ' "$log_dir/trace.tsv")" || fail "could not parse quasi task keys from the trace"

                [[ "$trace_matches" -eq 1 ]] ||
                    fail "expected one completed task for $expected_tag, found $trace_matches"
            done
        done
    done
}

# Verify the opt-in SVC run without adding any work to the default non-ML smoke path.
verify_svc_trace_keys() {
    basic_matches="$(awk -F '\t' '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "name") name_column = i
                if ($i == "status") status_column = i
            }
            if (!name_column || !status_column) exit 2
            next
        }
        index($(name_column), "SVC_BASIC") && $(status_column) == "COMPLETED" {
            matches += 1
        }
        END { print matches + 0 }
    ' "$log_dir/trace.tsv")" || fail "could not parse the SVC basic task from the trace"
    [[ "$basic_matches" -eq 1 ]] ||
        fail "expected one completed SVC_BASIC task, found $basic_matches"

    for ((replicate_id = 1; replicate_id <= smoke_nperms; replicate_id++)); do
        expected_tag="SVC_PERMUTE (replicate=$replicate_id)"
        permute_matches="$(awk -F '\t' -v expected="$expected_tag" '
            NR == 1 {
                for (i = 1; i <= NF; i++) {
                    if ($i == "name") name_column = i
                    if ($i == "status") status_column = i
                }
                if (!name_column || !status_column) exit 2
                next
            }
            index($(name_column), expected) && $(status_column) == "COMPLETED" {
                matches += 1
            }
            END { print matches + 0 }
        ' "$log_dir/trace.tsv")" || fail "could not parse the SVC permutation task from the trace"
        [[ "$permute_matches" -eq 1 ]] ||
            fail "expected one completed task for $expected_tag, found $permute_matches"
    done
}

verify_svc_outputs() {
    svc_output_dir="$output_dir/svctables"
    [[ -d "$svc_output_dir" ]] || fail "enabled SVC output directory is missing"
    [[ -f "$repo_dir/bin/zuzu_svc.py" ]] || fail "the SVC model helper is missing"

    svc_validation_dir="$smoke_root/svc-validation"
    svc_validation_artifacts="$svc_validation_dir/artifacts"
    svc_validation_sheets="$svc_validation_dir/sheets"
    mkdir -p "$svc_validation_artifacts" "$svc_validation_sheets"
    cp "$smoke_counts" "$svc_validation_dir/raw_counts.tsv"
    cp "$clean_counts" "$svc_validation_dir/analysis_counts.csv"
    cp "$clean_sheet" "$svc_validation_sheets/svc_basic.csv"

    svc_suffixes=(
        table.csv
        curve.csv
        folds.csv
        fs.png
        model.pickle
        metadata.json
        n_selected_features.txt
    )
    svc_prefixes=(svc_basic)
    for ((replicate_id = 1; replicate_id <= smoke_nperms; replicate_id++)); do
        svc_prefixes+=("svc_permute_replicate_$replicate_id")
    done

    for svc_prefix in "${svc_prefixes[@]}"; do
        for svc_suffix in "${svc_suffixes[@]}"; do
            svc_artifact="$svc_output_dir/${svc_prefix}_${svc_suffix}"
            require_file "$svc_artifact"
            cp -L "$svc_artifact" "$svc_validation_artifacts/"
        done
        require_png "$svc_output_dir/${svc_prefix}_fs.png"
    done

    staged_permutation_sheets=0
    while IFS= read -r metadata_file; do
        permutation_dir="${metadata_file%/*}"
        replicate_id="$(awk -F',' '
            NR == 1 {
                for (i = 1; i <= NF; i++) {
                    gsub(/^"|"$/, "", $i)
                    if ($i == "replicate") replicate_column = i
                }
                if (!replicate_column) exit 2
                next
            }
            NR == 2 {
                value = $(replicate_column)
                gsub(/^"|"$/, "", value)
                print value
                found = 1
                next
            }
            END { if (!found || NR != 2) exit 3 }
        ' "$metadata_file")" || fail "could not read SVC permutation replicate metadata"
        [[ "$replicate_id" =~ ^[1-9][0-9]*$ ]] ||
            fail "invalid SVC permutation replicate: $replicate_id"
        require_file "$permutation_dir/samplesheet_perm.csv"
        cp "$permutation_dir/samplesheet_perm.csv" \
            "$svc_validation_sheets/svc_permute_replicate_$replicate_id.csv"
        ((staged_permutation_sheets += 1))
    done < <(find "$work_dir" -mindepth 3 -maxdepth 3 -type f -name 'permutation_metadata.csv' -print | sort)
    [[ "$staged_permutation_sheets" -eq "$smoke_nperms" ]] ||
        fail "expected $smoke_nperms SVC permutation sample sheets, found $staged_permutation_sheets"

    svc_container="$(sed -n "s/^process.container[[:space:]]*=[[:space:]]*'\([^']*\)'[[:space:]]*$/\1/p" \
        "$repo_dir/nextflow.config" | head -n 1)"
    [[ -n "$svc_container" ]] || fail "could not determine the SVC container from nextflow.config"

    if ! docker run --rm --network none --read-only --tmpfs /tmp \
        -e PYTHONDONTWRITEBYTECODE=1 \
        -v "$svc_validation_dir:/validation:ro" \
        -v "$repo_dir/bin:/zuzu-bin:ro" \
        -v "$repo_dir/tests:/zuzu-tests:ro" \
        "$svc_container" \
        python3 /zuzu-tests/verify_svc_smoke.py \
        "$smoke_nperms" "$smoke_reflevel" "$smoke_altlevel"
    then
        fail "SVC artifact or saved-model validation failed"
    fi
}

for smoke_pdf in \
    QC_plots.pdf \
    permute_plot.pdf \
    upSet_basic.pdf \
    basic_overlap_violins.pdf \
    permute_hist.pdf \
    quasi_plot.pdf
do
    require_pdf "$output_dir/$smoke_pdf"
done

require_file "$output_dir/deseq_poorfits.RData"
require_file "$output_dir/quasi_results.csv"

# KI-027: with fewer than five permutations this is currently RData carrying a
# .pdf extension. Keep the default smoke run fast without blessing that format.
if [[ "$smoke_nperms" -ge 5 ]]; then
    require_pdf "$output_dir/permute_violins.pdf"
else
    require_file "$output_dir/permute_violins.pdf"
fi

if [[ "$smoke_synth" == 1 ]]; then
    require_pdf "$output_dir/synth_plot.pdf"
    require_file "$output_dir/synth_results.csv"
    require_file "$output_dir/synth_summary.csv"
else
    for disabled_synth_output in synth_plot.pdf synth_results.csv synth_summary.csv; do
        [[ ! -e "$output_dir/$disabled_synth_output" ]] ||
            fail "disabled synthetic output was published: $disabled_synth_output"
    done
fi
if [[ "$smoke_svc" == 0 ]]; then
    [[ ! -d "$output_dir/svctables" ]] || fail "disabled SVC output was published"
fi
require_file "$log_dir/trace.tsv"
verify_filter_warning
locate_clean_outputs
verify_common_filter
verify_reference_direction
verify_null_permutations
verify_permutation_trace_keys
verify_quasi_trace_keys
verify_quasi_anchor_contract
verify_quasi_metrics
if [[ "$smoke_svc" == 1 ]]; then
    verify_svc_trace_keys
    verify_svc_outputs
fi
if [[ "$smoke_synth" == 1 ]]; then
    verify_fullsynth_datasets
    verify_fullsynth_metrics
fi

# Materialize the valid plots before the temporary work tree and its published
# symlinks are removed. This directory is ignored by Git and refreshed in place.
mkdir -p "$plot_dir"
for smoke_pdf in \
    QC_plots.pdf \
    permute_plot.pdf \
    upSet_basic.pdf \
    basic_overlap_violins.pdf \
    permute_hist.pdf \
    quasi_plot.pdf
do
    cp -L "$output_dir/$smoke_pdf" "$plot_dir/$smoke_pdf"
done

if [[ "$smoke_synth" == 1 ]]; then
    cp -L "$output_dir/synth_plot.pdf" "$plot_dir/synth_plot.pdf"
else
    rm -f "$plot_dir/synth_plot.pdf"
fi

cp -L "$output_dir/quasi_results.csv" "$plot_dir/quasi_results.csv"

if [[ "$smoke_nperms" -ge 5 ]]; then
    cp -L "$output_dir/permute_violins.pdf" "$plot_dir/permute_violins.pdf"
else
    rm -f "$plot_dir/permute_violins.pdf"
fi

svc_plot_dir="$plot_dir/svc"
if [[ "$smoke_svc" == 1 ]]; then
    case "$svc_plot_dir" in
        "$plot_dir"/svc) rm -rf -- "$svc_plot_dir" ;;
        *) fail "refusing to refresh unexpected SVC artifact directory: $svc_plot_dir" ;;
    esac
    mkdir -p "$svc_plot_dir"
    for svc_artifact in "$output_dir"/svctables/svc_*; do
        [[ -e "$svc_artifact" ]] || fail "no SVC artifacts were available to retain"
        cp -L "$svc_artifact" "$svc_plot_dir/"
    done
else
    case "$svc_plot_dir" in
        "$plot_dir"/svc) rm -rf -- "$svc_plot_dir" ;;
        *) fail "refusing to remove unexpected SVC artifact directory: $svc_plot_dir" ;;
    esac
fi

printf 'Smoke-test artifacts retained at: %s\n' "$plot_dir"
