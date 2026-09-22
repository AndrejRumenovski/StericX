//! Corrected failure witnesses and descriptor-contract tests. The audit is read-only.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

const EXE: &str = env!("CARGO_BIN_EXE_stericx");

struct Workspace(PathBuf);

impl Workspace {
    fn new() -> Self {
        static SERIAL: AtomicUsize = AtomicUsize::new(0);
        let root = std::env::temp_dir().join(format!(
            "stericx_science_fix_{}_{}",
            std::process::id(),
            SERIAL.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&root).unwrap();
        Self(root)
    }

    fn file(&self, name: &str, content: &str) -> PathBuf {
        let path = self.0.join(name);
        std::fs::write(&path, content).unwrap();
        path
    }

    fn model(&self, aggregation: Option<&str>) -> PathBuf {
        let mut model: serde_json::Value =
            serde_json::from_str(include_str!("data/golden_training_report.json")).unwrap();
        if let Some(value) = aggregation {
            model["descriptor_aggregation"] = serde_json::json!(value);
        }
        self.file("model.json", &model.to_string())
    }

    fn reaction_csv(&self, weights: &str) -> PathBuf {
        self.file(
            "one.xyz",
            "3\nfirst conformer\nP 0 0 0\nC 0 0 1.8\nF 1.2 0 2.8\n",
        );
        self.file(
            "two.xyz",
            "3\nsecond conformer\nP 0 0 0\nC 0 0 1.8\nF 4.2 0 3.8\n",
        );
        self.file("reactions.csv", &format!(
            "Reaction_ID,Ligand_XYZ_Path,Attach_Atom_Idx,Primary_Bond_Vector_Idx,NBO_Charge,IR_Frequency,Temp_K,Exp_ddG_kcal_mol,Conformer_XYZ_Paths,Conformer_Boltzmann_Weights\nR,one.xyz,0,1,-0.35,1650,353.15,1.2,one.xyz;two.xyz,{weights}\n"
        ))
    }
}

impl Drop for Workspace {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn screen(model: &Path, library: &Path) -> Output {
    Command::new(EXE)
        .arg("screen")
        .arg(model)
        .arg("--library")
        .arg(library)
        .args(["--ascending", "--format", "json"])
        .output()
        .unwrap()
}

#[test]
fn original_nan_infinity_and_mismatch_evaluation_witnesses_are_rejected() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/model_evaluation");
    let temp = Workspace::new();
    for kind in ["nan", "infinity", "finite_mismatch"] {
        let predictions = root.join(format!("evaluate_{kind}_predictions.csv"));
        assert!(
            predictions.is_file(),
            "original witness missing: {}",
            predictions.display()
        );
        let output_path = temp.0.join(format!("{kind}.json"));
        let output = Command::new(EXE)
            .arg("evaluate")
            .arg("--data")
            .arg(root.join("synthetic_linear.sigpack"))
            .arg("--metadata")
            .arg(root.join("synthetic_linear_labels.csv"))
            .arg("--model")
            .arg(root.join("report.json"))
            .arg("--predictions")
            .arg(predictions)
            .arg("--output")
            .arg(&output_path)
            .output()
            .unwrap();
        assert!(!output.status.success(), "{kind} was accepted");
        assert!(
            !output_path.exists(),
            "invalid evaluation wrote successful metrics"
        );
        let expected = if kind == "finite_mismatch" {
            "does not match the supplied model"
        } else {
            "must both be finite"
        };
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains(expected), "{kind}: {stderr}");
    }
}

#[test]
fn weighted_screen_uses_exactly_the_descriptors_written_by_parse() {
    let temp = Workspace::new();
    let csv = temp.reaction_csv("1;3");
    let packed = temp.0.join("records.sigpack");
    let parse = Command::new(EXE)
        .arg("parse")
        .arg("--csv")
        .arg(&csv)
        .arg("--xyz-dir")
        .arg(&temp.0)
        .arg("--output")
        .arg(&packed)
        .output()
        .unwrap();
    assert!(
        parse.status.success(),
        "{}",
        String::from_utf8_lossy(&parse.stderr)
    );
    let model = temp.model(Some("supplied_weight_mean"));
    let output = screen(&model, &csv);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let matrix = steric_x::SigPackReader::open(&packed).unwrap();
    let record = matrix.records()[0];
    for descriptor in report["hits"][0]["descriptors"].as_array().unwrap() {
        let expected = match descriptor["name"].as_str().unwrap() {
            "sterimol_l" => record.l,
            "sterimol_b1" => record.b1,
            "sterimol_b5" => record.b5,
            "nbo_charge" => record.nbo_charge,
            "ir_frequency" => record.ir_freq,
            name => panic!("unexpected descriptor {name}"),
        };
        assert_eq!(descriptor["value"].as_f64().unwrap(), f64::from(expected));
    }
    assert_eq!(report["descriptor_aggregation"], "supplied_weight_mean");
    assert_ne!(report["hits"][0]["trust"], "reliable");
}

#[test]
fn single_geometry_honors_explicit_row_axes_and_identifies_supplied_mixtures() {
    let temp = Workspace::new();
    let csv = temp.reaction_csv("1");
    let source = std::fs::read_to_string(&csv)
        .unwrap()
        .replace("one.xyz;two.xyz", "one.xyz")
        .replace("R,one.xyz,0,1,", "R,one.xyz,1,0,");
    std::fs::write(&csv, &source).unwrap();
    let output = screen(&temp.model(Some("single_geometry")), &csv);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let molecule = steric_x::parse_coordinate_file(temp.0.join("one.xyz")).unwrap();
    let expected = steric_x::SterimolCalculator::compute(&molecule[0], 1, 0).unwrap();
    for descriptor in report["hits"][0]["descriptors"].as_array().unwrap() {
        let value = match descriptor["name"].as_str().unwrap() {
            "sterimol_l" => expected.l,
            "sterimol_b1" => expected.b1,
            "sterimol_b5" => expected.b5,
            _ => continue,
        };
        assert_eq!(descriptor["value"].as_f64().unwrap(), f64::from(value));
        assert_eq!(
            descriptor["source"],
            "computed_from_geometry_under_declared_aggregation"
        );
    }
    // A partial axis must not fall through to automatic donor selection.
    std::fs::write(&csv, source.replace("R,one.xyz,1,0,", "R,one.xyz,1,,")).unwrap();
    let invalid = screen(&temp.model(Some("single_geometry")), &csv);
    assert!(!invalid.status.success());
    assert!(String::from_utf8_lossy(&invalid.stderr).contains("requires both"));
    let mut lines = source.lines();
    std::fs::write(
        &csv,
        format!(
            "{},sterimol_b5\n{},3.1\n",
            lines.next().unwrap(),
            lines.next().unwrap()
        ),
    )
    .unwrap();
    let mixed = screen(&temp.model(Some("single_geometry")), &csv);
    assert!(mixed.status.success());
    let report: serde_json::Value = serde_json::from_slice(&mixed.stdout).unwrap();
    let descriptor = report["hits"][0]["descriptors"]
        .as_array()
        .unwrap()
        .iter()
        .find(|descriptor| descriptor["name"] == "sterimol_b5")
        .unwrap();
    assert_eq!(descriptor["value"].as_f64().unwrap(), f64::from(3.1_f32));
    assert_eq!(descriptor["source"], "supplied_value_caller_provenance");
}

#[test]
fn geometry_substitution_requires_the_declared_aggregation_and_explicit_weights() {
    for aggregation in [
        None,
        Some("supplied_record_values"),
        Some("single_geometry"),
    ] {
        let temp = Workspace::new();
        let csv = temp.reaction_csv("1;3");
        let output = screen(&temp.model(aggregation), &csv);
        assert!(
            !output.status.success(),
            "mismatched method {aggregation:?} was accepted"
        );
        let error = String::from_utf8_lossy(&output.stderr);
        assert!(
            error.contains("precomputed") || error.contains("multiple conformers"),
            "{error}"
        );
    }
    for weights in ["", "0;0", "-1;2", "NaN;2", "1"] {
        let temp = Workspace::new();
        let csv = temp.reaction_csv(weights);
        let output = screen(&temp.model(Some("supplied_weight_mean")), &csv);
        assert!(
            !output.status.success(),
            "invalid weights {weights:?} were accepted"
        );
    }
}
