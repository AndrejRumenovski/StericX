//! The response annotation must not invent stereochemical labels from numeric targets.

use std::path::PathBuf;
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

const EXE: &str = env!("CARGO_BIN_EXE_stericx");
const MAGNITUDE: &str = "Magnitude |ddG| from ddG_abs; larger values mean greater enantioselectivity; no R/S assignment.";

struct FitRun(PathBuf);

impl FitRun {
    fn new() -> Self {
        static COUNTER: AtomicUsize = AtomicUsize::new(0);
        let path = std::env::temp_dir().join(format!(
            "stericx_response_{}_{}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }

    fn fit(&self, convention: Option<&str>) -> Output {
        let mut command = Command::new(EXE);
        command.current_dir(env!("CARGO_MANIFEST_DIR")).args([
            "fit",
            "--data",
            "data/reactions.sigpack",
            "--metadata",
            "data/reactions_raw.csv",
            "--bootstrap",
            "20",
            "--permutations",
            "20",
            "--seed",
            "20260916",
            "--optimize",
            "maximize",
        ]);
        for (flag, filename) in [
            ("--output", "report.json"),
            ("--predictions", "predictions.csv"),
            ("--portable-model", "model.json"),
        ] {
            command.arg(flag).arg(self.0.join(filename));
        }
        if let Some(value) = convention {
            command.args(["--response-sign-convention", value]);
        }
        command.output().unwrap()
    }

    fn json(&self, name: &str) -> serde_json::Value {
        serde_json::from_slice(&std::fs::read(self.0.join(name)).unwrap()).unwrap()
    }
}

impl Drop for FitRun {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[test]
fn explicit_magnitude_is_inspectable_without_changing_fit_or_predictions() {
    let unspecified = FitRun::new();
    let magnitude = FitRun::new();
    for output in [unspecified.fit(None), magnitude.fit(Some(MAGNITUDE))] {
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
    let default_model = unspecified.json("model.json");
    assert!(
        default_model["inference"]["response"]["sign_convention"]
            .as_str()
            .unwrap()
            .starts_with("Not specified;")
    );
    assert_eq!(
        magnitude.json("model.json")["inference"]["response"]["sign_convention"],
        MAGNITUDE
    );
    assert_eq!(
        unspecified.json("report.json"),
        magnitude.json("report.json")
    );
    assert_eq!(
        std::fs::read(unspecified.0.join("predictions.csv")).unwrap(),
        std::fs::read(magnitude.0.join("predictions.csv")).unwrap()
    );
    let inspected = Command::new(EXE)
        .args(["model", "inspect"])
        .arg(magnitude.0.join("model.json"))
        .output()
        .unwrap();
    assert!(inspected.status.success());
    assert!(String::from_utf8_lossy(&inspected.stdout).contains(MAGNITUDE));
}

#[test]
fn empty_convention_fails_before_any_artifacts_are_written() {
    let run = FitRun::new();
    let output = run.fit(Some("   "));
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("target definition"));
    assert_eq!(std::fs::read_dir(&run.0).unwrap().count(), 0);
}
