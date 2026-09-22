//! Parallel file scheduling must preserve the serial CLI's observable contract.
//!
//! These tests reuse single-file outputs without introducing new scientific
//! expectations. The separately frozen CLI oracle anchors the numerical values.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

const EMPTY_BATCH_ERROR: &str = "error: no ligand files could be featurized\n";

struct Fixtures {
    root: PathBuf,
    large: PathBuf,
    small: PathBuf,
    malformed: PathBuf,
    no_donor: PathBuf,
    missing: PathBuf,
}

impl Fixtures {
    fn new() -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let root = std::env::temp_dir().join(format!(
            "stericx_descriptor_batch_{}_{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&root).unwrap();
        let fixtures = Self {
            large: root.join("large, ligand.xyz"),
            small: root.join("small.xyz"),
            malformed: root.join("truncated.xyz"),
            no_donor: root.join("no_donor.xyz"),
            missing: root.join("missing.xyz"),
            root,
        };
        let repository = Path::new(env!("CARGO_MANIFEST_DIR"));
        fs::copy(
            repository.join("data/xyz/SIG-NIHDA-785_fe24eb0f.xyz"),
            &fixtures.large,
        )
        .unwrap();
        fs::copy(
            repository.join("data/xyz/SIG-NIHDA-401_9d42bff1.xyz"),
            &fixtures.small,
        )
        .unwrap();
        fs::write(&fixtures.malformed, "2\ntruncated\nP 0 0 0\n").unwrap();
        fs::write(&fixtures.no_donor, "1\nno donor\nC 0 0 0\n").unwrap();
        fixtures
    }
}

impl Drop for Fixtures {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.root).unwrap();
    }
}

fn run(inputs: &[&Path], format: &str, threads: usize, extra: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_stericx"))
        .arg("descriptors")
        .args(inputs)
        .args(["--format", format])
        .args(extra)
        .env("RAYON_NUM_THREADS", threads.to_string())
        .env("LC_ALL", "C")
        .env_remove("STERICX_PROFILE_PATH")
        .output()
        .expect("descriptor CLI runs")
}

fn single_success(path: &Path, format: &str) -> String {
    let output = run(&[path], format, 1, &[]);
    assert!(output.status.success(), "{:?}", output);
    assert!(output.stderr.is_empty());
    String::from_utf8(output.stdout).unwrap()
}

fn single_skip(path: &Path) -> String {
    let output = run(&[path], "json", 1, &[]);
    assert_eq!(output.status.code(), Some(2));
    assert!(output.stdout.is_empty());
    let stderr = String::from_utf8(output.stderr).unwrap();
    let skip = stderr.strip_suffix(EMPTY_BATCH_ERROR).unwrap();
    assert!(skip.starts_with(&format!("skipped {}: ", path.display())));
    skip.to_owned()
}

fn combined_stdout(outputs: &[&str], format: &str) -> String {
    match format {
        "json" => format!(
            "[\n{}\n]\n",
            outputs
                .iter()
                .map(|output| output
                    .strip_prefix("[\n")
                    .unwrap()
                    .strip_suffix("\n]\n")
                    .unwrap())
                .collect::<Vec<_>>()
                .join(",\n")
        ),
        "text" => outputs.join("\n"),
        "csv" => {
            let header = outputs[0].split_once('\n').unwrap().0;
            let mut combined = format!("{header}\n");
            for output in outputs {
                let (actual_header, row) = output.split_once('\n').unwrap();
                assert_eq!(actual_header, header);
                combined.push_str(row);
            }
            combined
        }
        _ => unreachable!(),
    }
}

#[test]
fn mixed_batch_preserves_values_input_order_duplicates_and_skip_order_in_all_formats() {
    let f = Fixtures::new();
    let malformed = single_skip(&f.malformed);
    let no_donor = single_skip(&f.no_donor);
    let expected_stderr =
        format!("{malformed}{no_donor}{no_donor}featurized 4 of 7 files (3 skipped)\n");
    let inputs: [&Path; 7] = [
        &f.malformed,
        &f.large,
        &f.no_donor,
        &f.small,
        &f.large,
        &f.no_donor,
        &f.small,
    ];
    for format in ["json", "text", "csv"] {
        let large = single_success(&f.large, format);
        let small = single_success(&f.small, format);
        let expected_stdout = combined_stdout(&[&large, &small, &large, &small], format);
        for threads in [1, 2, 4, 6] {
            let output = run(&inputs, format, threads, &[]);
            assert_eq!(output.status.code(), Some(0), "{format}, {threads} threads");
            assert_eq!(
                output.stdout,
                expected_stdout.as_bytes(),
                "{format}, {threads} threads"
            );
            assert_eq!(
                output.stderr,
                expected_stderr.as_bytes(),
                "{format}, {threads} threads"
            );
        }
    }
}

#[test]
fn all_failed_inputs_emit_every_skip_in_order_then_one_terminal_error() {
    let f = Fixtures::new();
    let malformed = single_skip(&f.malformed);
    let no_donor = single_skip(&f.no_donor);
    let missing = single_skip(&f.missing);
    let expected_stderr = format!("{malformed}{malformed}{no_donor}{missing}{EMPTY_BATCH_ERROR}");
    for format in ["json", "text", "csv"] {
        for threads in [1, 2, 4, 6] {
            let output = run(
                &[&f.malformed, &f.malformed, &f.no_donor, &f.missing],
                format,
                threads,
                &[],
            );
            assert_eq!(output.status.code(), Some(2));
            assert!(output.stdout.is_empty(), "{format}, {threads} threads");
            assert_eq!(
                output.stderr,
                expected_stderr.as_bytes(),
                "{format}, {threads} threads"
            );
        }
    }
}

#[test]
fn batch_explicit_indices_are_rejected_before_reading_any_input() {
    let f = Fixtures::new();
    for format in ["json", "text", "csv"] {
        for threads in [1, 2, 4, 6] {
            for flag in ["--donor-index", "--reference-index"] {
                let output = run(
                    &[&f.missing, &f.malformed, &f.large],
                    format,
                    threads,
                    &[flag, "0"],
                );
                assert_eq!(output.status.code(), Some(2));
                assert!(output.stdout.is_empty());
                assert_eq!(
                output.stderr,
                b"error: --donor-index/--reference-index apply to a single file; omit it for batch runs\n"
            );
            }
        }
    }
}

#[test]
fn mixed_ensemble_batch_keeps_checked_means_and_single_file_results() {
    let f = Fixtures::new();
    let ensemble = f.root.join("finite ensemble.sdf");
    fs::copy(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("docs/scientific_remediation/geometry/finite_mean/finite_mean_4.sdf"),
        &ensemble,
    )
    .unwrap();
    let extra = ["--sterimol-axis", "coordination", "--density", "1"];
    for format in ["json", "text", "csv"] {
        let single = run(&[&ensemble], format, 1, &extra);
        let small = run(&[&f.small], format, 1, &extra);
        assert!(single.status.success() && single.stderr.is_empty());
        assert!(small.status.success() && small.stderr.is_empty());
        let expected = combined_stdout(
            &[
                std::str::from_utf8(&single.stdout).unwrap(),
                std::str::from_utf8(&small.stdout).unwrap(),
                std::str::from_utf8(&single.stdout).unwrap(),
            ],
            format,
        );
        for threads in [1, 2, 4, 6] {
            let again = run(&[&ensemble], format, threads, &extra);
            assert_eq!(again.status.code(), single.status.code());
            assert_eq!(again.stdout, single.stdout);
            assert_eq!(again.stderr, single.stderr);
            let batch = run(&[&ensemble, &f.small, &ensemble], format, threads, &extra);
            assert!(batch.status.success());
            assert!(batch.stderr.is_empty());
            assert_eq!(batch.stdout, expected.as_bytes());
        }
    }
}
