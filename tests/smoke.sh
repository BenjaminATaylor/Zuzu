#!/usr/bin/env bash

set -Eeuo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo_dir="$(cd "$script_dir/.." && pwd -P)"
source_sheet="$repo_dir/input/synthtest/synthsheet.csv"
source_counts="$repo_dir/input/synthtest/synthcounts.csv"
plot_dir="$script_dir/.smoke-plots"
smoke_nperms="${SMOKE_NPERMS:-1}"

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
# KI-003 regression detects accidental reuse or loss of either group size.
# CLEANINPUTS must leave both source files byte-for-byte unchanged.
smoke_sheet="$input_dir/samplesheet.csv"
smoke_counts="$input_dir/countsframe_clean.csv"
head -n 55 "$source_sheet" > "$smoke_sheet"

# The compact tracked fixture is CSV, while CLEANINPUTS expects tab-delimited counts.
# It has no embedded commas; keep the gene column plus samples 1-54 and convert
# the result mechanically to tab-delimited input.
awk -F',' 'BEGIN { OFS = "\t" }
    {
        for (i = 1; i <= 55; i++) {
            printf "%s%s", $i, (i == 55 ? ORS : OFS)
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
        --reflevel foo \
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
