#![cfg(feature = "profiling")]

use std::path::Path;
use std::process::Command;

#[test]
fn diagnostic_report_partitions_serial_wall_time_and_preserves_stdout() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let report =
        std::env::temp_dir().join(format!("stericx_profile_cli_{}.json", std::process::id()));
    let run = |enabled: bool| {
        let mut command = Command::new(env!("CARGO_BIN_EXE_stericx"));
        command
            .current_dir(root)
            .args([
                "descriptors",
                "data/xyz/SIG-NIHDA-401_9d42bff1.xyz",
                "--format",
                "json",
            ])
            .env_remove("STERICX_PROFILE_PATH");
        if enabled {
            command.env("STERICX_PROFILE_PATH", &report);
        }
        command.output().unwrap()
    };
    let plain = run(false);
    let profiled = run(true);
    assert!(plain.status.success());
    assert!(profiled.status.success());
    assert_eq!(plain.stdout, profiled.stdout);
    assert_eq!(plain.stderr, profiled.stderr);
    let data: serde_json::Value = serde_json::from_slice(&std::fs::read(&report).unwrap()).unwrap();
    std::fs::remove_file(report).unwrap();
    assert_eq!(data["dropped_scopes"], 0);
    assert_eq!(data["observed_threads"], 1);
    let functions = data["functions"].as_array().unwrap();
    let root = functions
        .iter()
        .find(|entry| entry["function"] == "cli::run")
        .unwrap();
    let root_ns = root["inclusive_ns"].as_u64().unwrap();
    let exclusive_ns: u64 = functions
        .iter()
        .map(|entry| entry["exclusive_ns"].as_u64().unwrap())
        .sum();
    assert_eq!(exclusive_ns, root_ns);
    assert!(root_ns <= data["session_wall_ns"].as_u64().unwrap());
    for stage in [
        "file_parsing",
        "donor_bond_detection",
        "sterimol",
        "buried_volume",
        "pyramidalization",
        "conformer_processing",
        "output",
    ] {
        assert!(
            data["stage_exclusive_ns"][stage].as_u64().unwrap() > 0,
            "{stage}"
        );
    }
    assert!(data["allocations"]["allocation_calls"].as_u64().unwrap() > 0);
    assert!(
        data["allocations"]["peak_live_requested_bytes"]
            .as_u64()
            .unwrap()
            > 0
    );
}
