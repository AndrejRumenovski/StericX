//! Generic output helpers shared by the command handlers: atomic file writes,
//! resident-memory and timing metrics, and small formatting utilities.

use serde::Serialize;
use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Duration;

fn temporary_sibling(path: &Path) -> PathBuf {
    let mut value = path.as_os_str().to_os_string();
    value.push(".tmp");
    PathBuf::from(value)
}

pub(crate) fn ensure_parent(path: &Path) -> std::io::Result<()> {
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)?;
    }
    Ok(())
}

pub(crate) fn atomic_write_json<T: Serialize>(
    value: &T,
    path: &Path,
) -> Result<(), Box<dyn Error>> {
    ensure_parent(path)?;
    let temporary = temporary_sibling(path);
    let json = serde_json::to_string_pretty(value)?;
    fs::write(&temporary, format!("{json}\n"))?;
    fs::rename(temporary, path)?;
    Ok(())
}

pub(crate) fn atomic_write_csv_rows<T: Serialize>(
    rows: &[T],
    path: &Path,
) -> Result<(), Box<dyn Error>> {
    ensure_parent(path)?;
    let temporary = temporary_sibling(path);
    {
        let mut writer = csv::Writer::from_path(&temporary)?;
        for row in rows {
            writer.serialize(row)?;
        }
        writer.flush()?;
    }
    fs::rename(temporary, path)?;
    Ok(())
}

pub(crate) fn display_optional_metric(value: Option<f64>) -> String {
    value.map_or_else(|| "unavailable".into(), |value| format!("{value:.8}"))
}

pub(crate) fn fnv1a64(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |hash, byte| {
        (hash ^ u64::from(*byte)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

pub(crate) fn resident_memory_bytes() -> Option<u64> {
    let status = fs::read_to_string("/proc/self/status").ok()?;
    let line = status.lines().find(|line| line.starts_with("VmRSS:"))?;
    let kibibytes = line.split_whitespace().nth(1)?.parse::<u64>().ok()?;
    kibibytes.checked_mul(1_024)
}

pub(crate) fn print_memory_metrics(start: Option<u64>, end: Option<u64>) {
    match start {
        Some(bytes) => println!("rss_start_bytes={bytes}"),
        None => println!("rss_start_bytes=unavailable"),
    }
    match end {
        Some(bytes) => println!("rss_end_bytes={bytes}"),
        None => println!("rss_end_bytes=unavailable"),
    }
    match (start, end) {
        (Some(start), Some(end)) => {
            println!("rss_delta_bytes={}", i128::from(end) - i128::from(start));
        }
        _ => println!("rss_delta_bytes=unavailable"),
    }
}

pub(crate) fn millis(duration: Duration) -> f64 {
    duration.as_secs_f64() * 1_000.0
}
