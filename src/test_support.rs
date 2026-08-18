//! Shared helpers for the binary's unit tests.

use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

/// A unique, process- and time-stamped scratch directory path for a test.
pub(crate) fn temporary_directory(test_name: &str) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!(
        "steric_x_cli_{test_name}_{}_{nonce}",
        std::process::id()
    ))
}
