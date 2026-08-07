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

cleanup() {
    smoke_status=$?
    trap - EXIT

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

# These names deliberately avoid CLEANINPUTS' current fixed output basenames.
# Once KI-001 is fixed, use those collision-prone names here as a regression check.
smoke_sheet="$input_dir/smoke_sheet_input.csv"
smoke_counts="$input_dir/smoke_counts_input.tsv"
cp "$source_sheet" "$smoke_sheet"

# The compact tracked fixture is CSV, while CLEANINPUTS expects tab-delimited counts.
# This fixture has no embedded commas, so a mechanical conversion is sufficient.
tr ',' '\t' < "$source_counts" > "$smoke_counts"

source_sheet_before="$(cksum < "$source_sheet")"
source_counts_before="$(cksum < "$source_counts")"
smoke_sheet_before="$(cksum < "$smoke_sheet")"
smoke_counts_before="$(cksum < "$smoke_counts")"

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
        --truecutoff 1e-50 \
        --synthstep false \
        --mlstep false \
        --pareml true \
        -w "$work_dir" \
        -with-trace "$log_dir/trace.tsv"
) 2>&1 | tee "$log_dir/console.log"

[[ "$(cksum < "$source_sheet")" == "$source_sheet_before" ]] || fail "source sample sheet changed"
[[ "$(cksum < "$source_counts")" == "$source_counts_before" ]] || fail "source count matrix changed"
[[ "$(cksum < "$smoke_sheet")" == "$smoke_sheet_before" ]] || fail "staged sample sheet changed"
[[ "$(cksum < "$smoke_counts")" == "$smoke_counts_before" ]] || fail "staged count matrix changed"

require_file() {
    [[ -s "$1" ]] || fail "missing or empty artifact: $1"
}

require_pdf() {
    require_file "$1"
    LC_ALL=C head -c 5 "$1" | grep -q '^%PDF-' || fail "not a valid PDF: $1"
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
