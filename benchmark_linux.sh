#!/usr/bin/env bash
#
# Native Linux benchmark harness for StericX.
#
# Workloads:
#   1. Sterimol extraction from 10,000 XYZ files.
#   2. Flat binary packing of those 10,000 records.
#   3. Rayon/AVX2 regression over 1,000,000 memory-mapped records.
#
# GNU time observes the complete `parse` process, which contains both workload
# 1 and workload 2. StericX's own phase timers provide their separate latency;
# RSS, CPU, and page-fault figures are intentionally reported as shared process
# metrics for those two phases.

set -Eeuo pipefail
export LC_ALL=C

readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly BINARY="${SCRIPT_DIR}/target/release/stericx"
readonly SOURCE_CSV="${SCRIPT_DIR}/data/reactions_raw.csv"
readonly SOURCE_XYZ_DIR="${SCRIPT_DIR}/data/xyz"
readonly RESULTS_JSON="${SCRIPT_DIR}/docs/benchmark_results.json"
readonly XYZ_RECORD_COUNT=10000
readonly PREDICTION_RECORD_COUNT=1000000
readonly RECORD_BYTES=64
readonly THREAD_COUNT="${RAYON_NUM_THREADS:-$(nproc)}"

BENCH_ROOT=""

die() {
    printf 'benchmark error: %s\n' "$*" >&2
    exit 1
}

cleanup() {
    if [[ -n "${BENCH_ROOT}" &&
          -d "${BENCH_ROOT}" &&
          "${BENCH_ROOT}" == "${SCRIPT_DIR}"/target/stericx-benchmark.* ]]; then
        rm -rf -- "${BENCH_ROOT}"
    fi
}
trap cleanup EXIT INT TERM

require_command() {
    command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"
}

metric_from_gnu_time() {
    local metric="$1"
    local file="$2"
    awk -v metric="${metric}" '
        index($0, metric) {
            sub(/^.*: /, "")
            print
            exit
        }
    ' "${file}"
}

elapsed_to_seconds() {
    local elapsed="$1"
    awk -v elapsed="${elapsed}" 'BEGIN {
        count = split(elapsed, fields, ":")
        if (count == 3) {
            printf "%.6f", fields[1] * 3600 + fields[2] * 60 + fields[3]
        } else if (count == 2) {
            printf "%.6f", fields[1] * 60 + fields[2]
        } else {
            printf "%.6f", fields[1]
        }
    }'
}

milliseconds_to_seconds() {
    awk -v milliseconds="$1" 'BEGIN { printf "%.6f", milliseconds / 1000.0 }'
}

records_per_second() {
    awk -v records="$1" -v seconds="$2" 'BEGIN {
        if (seconds > 0) {
            printf "%.1f", records / seconds
        } else {
            printf "0.0"
        }
    }'
}

kib_to_mib() {
    awk -v kib="$1" 'BEGIN { printf "%.2f", kib / 1024.0 }'
}

internal_metric() {
    local key="$1"
    local file="$2"
    awk -F= -v key="${key}" '$1 == key { print $2; exit }' "${file}"
}

json_escape() {
    sed -e 's/\\/\\\\/g' -e 's/"/\\"/g' <<<"$1" | tr -d '\n'
}

run_timed() {
    local stdout_file="$1"
    local time_file="$2"
    shift 2
    if ! /usr/bin/time -v -o "${time_file}" "$@" >"${stdout_file}" 2>&1; then
        printf '\nFailed command output:\n' >&2
        cat "${stdout_file}" >&2
        printf '\nGNU time output:\n' >&2
        cat "${time_file}" >&2
        return 1
    fi
}

parse_time_metrics() {
    local file="$1"
    local prefix="$2"
    local elapsed cpu rss major minor
    elapsed="$(metric_from_gnu_time 'Elapsed (wall clock)' "${file}")"
    cpu="$(metric_from_gnu_time 'Percent of CPU this job got' "${file}")"
    rss="$(metric_from_gnu_time 'Maximum resident set size' "${file}")"
    major="$(metric_from_gnu_time 'Major (requiring I/O) page faults' "${file}")"
    minor="$(metric_from_gnu_time 'Minor (reclaiming a frame) page faults' "${file}")"
    [[ -n "${elapsed}" && -n "${cpu}" && -n "${rss}" &&
       -n "${major}" && -n "${minor}" ]] ||
        die "could not parse GNU time metrics from ${file}"

    printf -v "${prefix}_wall_seconds" '%s' "$(elapsed_to_seconds "${elapsed}")"
    printf -v "${prefix}_cpu_percent" '%s' "${cpu%\%}"
    printf -v "${prefix}_rss_kb" '%s' "${rss}"
    printf -v "${prefix}_major_faults" '%s' "${major}"
    printf -v "${prefix}_minor_faults" '%s' "${minor}"
}

generate_xyz_workload() {
    local destination_dir="$1"
    local destination_csv="$2"
    local -a source_paths=()
    local -a attach_indices=()
    local -a neighbor_indices=()
    local -a nbo_values=()
    local -a ir_values=()
    local -a temp_values=()
    local -a ddg_values=()
    local reaction_id xyz_path attach neighbor nbo ir temp ddg

    exec 3<"${SOURCE_CSV}"
    IFS= read -r _header <&3
    while IFS=, read -r reaction_id xyz_path attach neighbor nbo ir temp ddg _rest <&3; do
        ddg="${ddg%$'\r'}"
        [[ -n "${xyz_path}" ]] || continue
        local source_path="${SOURCE_XYZ_DIR}/$(basename -- "${xyz_path}")"
        [[ -f "${source_path}" ]] ||
            die "CSV references missing XYZ source: ${source_path}"
        source_paths+=("${source_path}")
        attach_indices+=("${attach}")
        neighbor_indices+=("${neighbor}")
        nbo_values+=("${nbo}")
        ir_values+=("${ir}")
        temp_values+=("${temp}")
        ddg_values+=("${ddg}")
    done
    exec 3<&-
    local source_count="${#source_paths[@]}"
    (( source_count > 0 )) || die "source CSV contains no usable rows"

    mkdir -p -- "${destination_dir}"
    exec 4>"${destination_csv}"
    printf '%s\n' \
        'Reaction_ID,Ligand_XYZ_Path,Attach_Atom_Idx,Primary_Bond_Vector_Idx,NBO_Charge,IR_Frequency,Temp_K,Exp_ddG_kcal_mol' \
        >&4

    local index source_index filename
    for ((index = 0; index < XYZ_RECORD_COUNT; index++)); do
        source_index=$((index % source_count))
        printf -v filename 'bench_%05d.xyz' "${index}"
        if ! ln -- "${source_paths[source_index]}" "${destination_dir}/${filename}" \
            2>/dev/null; then
            cp -- "${source_paths[source_index]}" "${destination_dir}/${filename}"
        fi
        printf 'BENCH-%05d,%s,%s,%s,%s,%s,%s,%s\n' \
            "${index}" \
            "${filename}" \
            "${attach_indices[source_index]}" \
            "${neighbor_indices[source_index]}" \
            "${nbo_values[source_index]}" \
            "${ir_values[source_index]}" \
            "${temp_values[source_index]}" \
            "${ddg_values[source_index]}" \
            >&4
    done
    exec 4>&-
}

generate_prediction_workload() {
    local ten_thousand_pack="$1"
    local million_pack="$2"
    local repeats=$((PREDICTION_RECORD_COUNT / XYZ_RECORD_COUNT))
    (( PREDICTION_RECORD_COUNT % XYZ_RECORD_COUNT == 0 )) ||
        die "prediction record count must be divisible by XYZ record count"

    : >"${million_pack}"
    local iteration
    for ((iteration = 0; iteration < repeats; iteration++)); do
        cat -- "${ten_thousand_pack}" >>"${million_pack}"
    done
    local expected_bytes=$((PREDICTION_RECORD_COUNT * RECORD_BYTES))
    local actual_bytes
    actual_bytes="$(stat -c '%s' "${million_pack}")"
    [[ "${actual_bytes}" -eq "${expected_bytes}" ]] ||
        die "prediction pack has ${actual_bytes} bytes; expected ${expected_bytes}"
}

print_summary() {
    printf '\nStericX Linux benchmark summary\n'
    printf '%-22s %10s %11s %11s %14s %9s %12s %9s %9s\n' \
        'Phase' 'Records' 'Phase(s)' 'Process(s)' 'Records/s' 'CPU(%)' \
        'Peak RSS MB' 'Maj PF' 'Min PF'
    printf '%-22s %10s %11s %11s %14s %9s %12s %9s %9s\n' \
        '----------------------' '----------' '-----------' '-----------' \
        '--------------' '---------' '------------' '---------' '---------'
    printf '%-22s %10d %11s %11s %14s %9s %12s %9s %9s\n' \
        'Sterimol extraction' "${XYZ_RECORD_COUNT}" \
        "${extraction_phase_seconds}" "${parse_wall_seconds}" \
        "${extraction_throughput}" "${parse_cpu_percent}" \
        "$(kib_to_mib "${parse_rss_kb}")" "${parse_major_faults}" \
        "${parse_minor_faults}"
    printf '%-22s %10d %11s %11s %14s %9s %12s %9s %9s\n' \
        'Binary packing' "${XYZ_RECORD_COUNT}" \
        "${packing_phase_seconds}" "${parse_wall_seconds}" \
        "${packing_throughput}" "${parse_cpu_percent}" \
        "$(kib_to_mib "${parse_rss_kb}")" "${parse_major_faults}" \
        "${parse_minor_faults}"
    printf '%-22s %10d %11s %11s %14s %9s %12s %9s %9s\n' \
        'RegressX prediction' "${PREDICTION_RECORD_COUNT}" \
        "${prediction_phase_seconds}" "${predict_wall_seconds}" \
        "${prediction_throughput}" "${predict_cpu_percent}" \
        "$(kib_to_mib "${predict_rss_kb}")" "${predict_major_faults}" \
        "${predict_minor_faults}"
    printf '\nResults JSON: %s\n' "${RESULTS_JSON}"
}

main() {
    require_command cargo
    require_command rustc
    require_command nproc
    require_command awk
    require_command sed
    require_command stat
    [[ -x /usr/bin/time ]] || die 'GNU /usr/bin/time is required'
    [[ -f "${SOURCE_CSV}" ]] || die "source CSV not found: ${SOURCE_CSV}"
    [[ -d "${SOURCE_XYZ_DIR}" ]] || die "source XYZ directory not found"

    printf 'Building StericX with native CPU optimizations...\n'
    (
        cd -- "${SCRIPT_DIR}"
        RUSTFLAGS="-C target-cpu=native" cargo build --release
    )
    [[ -x "${BINARY}" ]] || die "release binary was not produced: ${BINARY}"

    mkdir -p -- "${SCRIPT_DIR}/target"
    BENCH_ROOT="$(mktemp -d "${SCRIPT_DIR}/target/stericx-benchmark.XXXXXX")"
    local xyz_benchmark_dir="${BENCH_ROOT}/xyz"
    local benchmark_csv="${BENCH_ROOT}/reactions_10000.csv"
    local pack_10000="${BENCH_ROOT}/reactions_10000.sigpack"
    local pack_1000000="${BENCH_ROOT}/reactions_1000000.sigpack"
    local weights_json="${BENCH_ROOT}/weights.json"
    local parse_stdout="${BENCH_ROOT}/parse.stdout"
    local parse_time="${BENCH_ROOT}/parse.time"
    local predict_stdout="${BENCH_ROOT}/predict.stdout"
    local predict_time="${BENCH_ROOT}/predict.time"

    printf 'Generating 10,000-file XYZ workload...\n'
    generate_xyz_workload "${xyz_benchmark_dir}" "${benchmark_csv}"

    printf 'Benchmarking Sterimol extraction and binary packing...\n'
    run_timed "${parse_stdout}" "${parse_time}" \
        env RAYON_NUM_THREADS="${THREAD_COUNT}" \
        "${BINARY}" parse \
        --csv "${benchmark_csv}" \
        --xyz-dir "${xyz_benchmark_dir}" \
        --output "${pack_10000}"
    parse_time_metrics "${parse_time}" parse

    local extraction_ms packing_ms
    extraction_ms="$(internal_metric geometry_compute_ms "${parse_stdout}")"
    packing_ms="$(internal_metric binary_export_ms "${parse_stdout}")"
    [[ -n "${extraction_ms}" && -n "${packing_ms}" ]] ||
        die "StericX parse output did not contain phase timings"
    extraction_phase_seconds="$(milliseconds_to_seconds "${extraction_ms}")"
    packing_phase_seconds="$(milliseconds_to_seconds "${packing_ms}")"
    extraction_throughput="$(
        records_per_second "${XYZ_RECORD_COUNT}" "${extraction_phase_seconds}"
    )"
    packing_throughput="$(
        records_per_second "${XYZ_RECORD_COUNT}" "${packing_phase_seconds}"
    )"

    printf 'Generating 1,000,000-record mapped prediction workload...\n'
    generate_prediction_workload "${pack_10000}" "${pack_1000000}"
    printf '%s\n' \
        '[0.10, 0.20, -0.10, 0.30, 0.50, -0.20, 0.10, 0.001]' \
        >"${weights_json}"

    printf 'Benchmarking multi-threaded SIMD RegressX prediction...\n'
    run_timed "${predict_stdout}" "${predict_time}" \
        env RAYON_NUM_THREADS="${THREAD_COUNT}" \
        "${BINARY}" predict \
        --data "${pack_1000000}" \
        --weights "${weights_json}"
    parse_time_metrics "${predict_time}" predict

    local prediction_ms
    prediction_ms="$(internal_metric prediction_latency_ms "${predict_stdout}")"
    [[ -n "${prediction_ms}" ]] ||
        die "StericX prediction output did not contain inference latency"
    prediction_phase_seconds="$(milliseconds_to_seconds "${prediction_ms}")"
    prediction_throughput="$(
        records_per_second \
            "${PREDICTION_RECORD_COUNT}" \
            "${prediction_phase_seconds}"
    )"

    mkdir -p -- "$(dirname -- "${RESULTS_JSON}")"
    local generated_at rust_version
    generated_at="$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
    rust_version="$(json_escape "$(rustc --version)")"

    cat >"${RESULTS_JSON}" <<JSON
{
  "generated_at_utc": "${generated_at}",
  "build": {
    "rustc": "${rust_version}",
    "rustflags": "-C target-cpu=native"
  },
  "methodology": {
    "record_bytes": ${RECORD_BYTES},
    "time_utility": "/usr/bin/time -v",
    "parse_process_metrics_shared": true,
    "note": "Sterimol extraction and packing are internal phases of one parse process; phase latency is isolated by StericX timers while GNU time CPU, RSS, and page faults are shared."
  },
  "workloads": {
    "sterimol_extraction": {
      "records": ${XYZ_RECORD_COUNT},
      "phase_seconds": ${extraction_phase_seconds},
      "process_wall_seconds": ${parse_wall_seconds},
      "throughput_records_per_second": ${extraction_throughput},
      "cpu_utilization_percent": ${parse_cpu_percent},
      "peak_rss_kb": ${parse_rss_kb},
      "peak_rss_mb": $(kib_to_mib "${parse_rss_kb}"),
      "major_page_faults": ${parse_major_faults},
      "minor_page_faults": ${parse_minor_faults}
    },
    "binary_packing": {
      "records": ${XYZ_RECORD_COUNT},
      "output_bytes": $(stat -c '%s' "${pack_10000}"),
      "phase_seconds": ${packing_phase_seconds},
      "process_wall_seconds": ${parse_wall_seconds},
      "throughput_records_per_second": ${packing_throughput},
      "cpu_utilization_percent": ${parse_cpu_percent},
      "peak_rss_kb": ${parse_rss_kb},
      "peak_rss_mb": $(kib_to_mib "${parse_rss_kb}"),
      "major_page_faults": ${parse_major_faults},
      "minor_page_faults": ${parse_minor_faults}
    },
    "regressx_prediction": {
      "records": ${PREDICTION_RECORD_COUNT},
      "mapped_bytes": $(stat -c '%s' "${pack_1000000}"),
      "rayon_threads": ${THREAD_COUNT},
      "phase_seconds": ${prediction_phase_seconds},
      "process_wall_seconds": ${predict_wall_seconds},
      "throughput_records_per_second": ${prediction_throughput},
      "cpu_utilization_percent": ${predict_cpu_percent},
      "peak_rss_kb": ${predict_rss_kb},
      "peak_rss_mb": $(kib_to_mib "${predict_rss_kb}"),
      "major_page_faults": ${predict_major_faults},
      "minor_page_faults": ${predict_minor_faults}
    }
  }
}
JSON

    print_summary
}

main "$@"
