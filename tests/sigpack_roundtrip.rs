use std::time::{SystemTime, UNIX_EPOCH};
use steric_x::{PackedReactionRecord, SigPackReader, SigPackWriter};

#[test]
fn round_trips_one_thousand_records_with_identical_floats() {
    let records: Vec<PackedReactionRecord> = (0..1_000)
        .map(|index| {
            let value = index as f32;
            PackedReactionRecord {
                l: value * 0.01,
                b1: value * 0.02 + 1.0,
                b5: value * 0.03 + 2.0,
                nbo_charge: value * -0.001,
                ir_freq: 1_500.0 + value,
                temp_k: 273.15 + value * 0.01,
                exp_ddg: value.sin(),
                reserved: [value; 9],
            }
        })
        .collect();
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = std::env::temp_dir().join(format!(
        "steric_x_integration_{}_{nonce}.sigpack",
        std::process::id()
    ));

    SigPackWriter::export(&records, &path).unwrap();
    let reader = SigPackReader::open(&path).unwrap();

    assert_eq!(reader.len(), 1_000);
    assert_eq!(reader.records(), records.as_slice());
    for (expected, actual) in records.iter().zip(reader.records()) {
        assert_eq!(actual.l.to_bits(), expected.l.to_bits());
        assert_eq!(actual.b1.to_bits(), expected.b1.to_bits());
        assert_eq!(actual.b5.to_bits(), expected.b5.to_bits());
        assert_eq!(actual.nbo_charge.to_bits(), expected.nbo_charge.to_bits());
        assert_eq!(actual.ir_freq.to_bits(), expected.ir_freq.to_bits());
        assert_eq!(actual.temp_k.to_bits(), expected.temp_k.to_bits());
        assert_eq!(actual.exp_ddg.to_bits(), expected.exp_ddg.to_bits());
    }

    drop(reader);
    std::fs::remove_file(path).unwrap();
}
