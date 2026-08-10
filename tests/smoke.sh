#!/usr/bin/env bash

set -Eeuo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo_dir="$(cd "$script_dir/.." && pwd -P)"
source_sheet="$repo_dir/input/synthtest/synthsheet.csv"
source_counts="$repo_dir/input/synthtest/synthcounts.csv"
plot_dir="$script_dir/.smoke-plots"
smoke_nperms="${SMOKE_NPERMS:-1}"
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

# The relaxed test-only anchor threshold makes the compact fixture exercise the
# quasi recovery calculations rather than only their too-few-anchors fallback.
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
        --synthstep false \
        --mlstep false \
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

verify_quasi_sample_sizes() {
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

        size_values="$(awk -F',' '
            NR == 1 {
                for (i = 1; i <= NF; i++) {
                    gsub(/^"|"$/, "", $i)
                    column[$i] = i
                }
                if (!("refnum" in column) || !("altnum" in column) || !("totalnum" in column)) {
                    exit 2
                }
                next
            }
            NR == 2 {
                for (i = 1; i <= NF; i++) {
                    gsub(/^"|"$/, "", $i)
                }
                refnum = $(column["refnum"])
                altnum = $(column["altnum"])
                totalnum = $(column["totalnum"])
                if (refnum !~ /^[0-9]+$/ || altnum !~ /^[0-9]+$/ || totalnum !~ /^[0-9]+$/) {
                    exit 3
                }
                if (totalnum != refnum + altnum) {
                    exit 4
                }
                printf "%d %d %d\n", refnum, altnum, totalnum
                next
            }
            END {
                if (NR != 2) {
                    exit 5
                }
            }
        ' "$quasi_result")" || fail "invalid quasi result metadata: $quasi_result"

        read -r observed_ref observed_alt _ <<< "$size_values"
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

# KI-027: with fewer than five permutations this is currently RData carrying a
# .pdf extension. Keep the default smoke run fast without blessing that format.
if [[ "$smoke_nperms" -ge 5 ]]; then
    require_pdf "$output_dir/permute_violins.pdf"
else
    require_file "$output_dir/permute_violins.pdf"
fi

[[ ! -e "$output_dir/synth_plot.pdf" ]] || fail "disabled synthetic output was published"
[[ ! -d "$output_dir/svctables" ]] || fail "disabled SVC output was published"
require_file "$log_dir/trace.tsv"
verify_filter_warning
locate_clean_outputs
verify_common_filter
verify_reference_direction
verify_null_permutations
verify_permutation_trace_keys
verify_quasi_trace_keys
verify_quasi_sample_sizes

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

if [[ "$smoke_nperms" -ge 5 ]]; then
    cp -L "$output_dir/permute_violins.pdf" "$plot_dir/permute_violins.pdf"
else
    rm -f "$plot_dir/permute_violins.pdf"
fi

printf 'Smoke-test plots retained at: %s\n' "$plot_dir"
