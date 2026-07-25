use super::{PackedReactionRecord, PackedReactionRecordV2, SigPackHeaderV2};
use memmap2::{Mmap, MmapOptions};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

/// Flat binary exporter for `.sigpack` record matrices.
#[derive(Clone, Copy, Debug, Default)]
pub struct SigPackWriter;

impl SigPackWriter {
    /// Serializes records as one contiguous native-endian byte stream.
    ///
    /// `PackedReactionRecord` is `Pod`, so the cast cannot expose padding or
    /// invalid bit patterns and requires no per-record serialization.
    pub fn export(records: &[PackedReactionRecord], output_path: &Path) -> io::Result<()> {
        let file = File::create(output_path)?;
        let mut writer = BufWriter::new(file);
        let bytes: &[u8] = bytemuck::cast_slice(records);
        writer.write_all(bytes)?;
        writer.flush()
    }
}

/// Read-only zero-copy view over a flat `.sigpack` record matrix.
pub struct SigPackReader {
    map: Option<Mmap>,
}

/// Version-two exporter with a validated 64-byte schema header.
#[derive(Clone, Copy, Debug, Default)]
pub struct SigPackV2Writer;

impl SigPackV2Writer {
    pub fn export(records: &[PackedReactionRecordV2], output_path: &Path) -> io::Result<()> {
        let file = File::create(output_path)?;
        let mut writer = BufWriter::new(file);
        let header = SigPackHeaderV2::new(records.len());
        writer.write_all(bytemuck::bytes_of(&header))?;
        writer.write_all(bytemuck::cast_slice(records))?;
        writer.flush()
    }
}

/// Zero-copy reader for version-two descriptor matrices.
pub struct SigPackV2Reader {
    map: Mmap,
}

impl SigPackV2Reader {
    pub fn open(input_path: &Path) -> io::Result<Self> {
        let file = File::open(input_path)?;
        let byte_len = usize::try_from(file.metadata()?.len()).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "sigpack v2 file is too large to address",
            )
        })?;
        if byte_len < size_of::<SigPackHeaderV2>() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "sigpack v2 file is shorter than its header",
            ));
        }
        // SAFETY: the file is opened read-only and Mmap owns the mapping.
        let map = unsafe { MmapOptions::new().map(&file)? };
        let header =
            bytemuck::try_from_bytes::<SigPackHeaderV2>(&map[..size_of::<SigPackHeaderV2>()])
                .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error.to_string()))?;
        validate_v2_header(header, byte_len)?;
        Ok(Self { map })
    }

    #[must_use]
    pub fn header(&self) -> &SigPackHeaderV2 {
        bytemuck::from_bytes(&self.map[..size_of::<SigPackHeaderV2>()])
    }

    #[must_use]
    pub fn records(&self) -> &[PackedReactionRecordV2] {
        bytemuck::cast_slice(&self.map[size_of::<SigPackHeaderV2>()..])
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.records().len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.records().is_empty()
    }
}

fn validate_v2_header(header: &SigPackHeaderV2, byte_len: usize) -> io::Result<()> {
    if header.magic != SigPackHeaderV2::MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "sigpack v2 magic is invalid",
        ));
    }
    if header.schema_version != SigPackHeaderV2::SCHEMA_VERSION {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "unsupported sigpack schema version {}",
                header.schema_version
            ),
        ));
    }
    if header.endian_marker != SigPackHeaderV2::ENDIAN_MARKER {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "sigpack v2 byte order does not match this host",
        ));
    }
    if header.record_size as usize != size_of::<PackedReactionRecordV2>() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "sigpack v2 record size does not match this build",
        ));
    }
    let record_count = usize::try_from(header.record_count).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "sigpack v2 record count is too large",
        )
    })?;
    let expected = size_of::<SigPackHeaderV2>()
        .checked_add(
            record_count
                .checked_mul(size_of::<PackedReactionRecordV2>())
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "sigpack v2 size overflows")
                })?,
        )
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "sigpack v2 size overflows"))?;
    if expected != byte_len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("sigpack v2 declares {expected} bytes but file contains {byte_len}"),
        ));
    }
    Ok(())
}

impl SigPackReader {
    /// Opens and memory-maps a `.sigpack` file.
    ///
    /// A non-empty file must contain an integral number of 64-byte records.
    pub fn open(input_path: &Path) -> io::Result<SigPackReader> {
        let file = File::open(input_path)?;
        let byte_len = usize::try_from(file.metadata()?.len()).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "sigpack file is too large to address",
            )
        })?;
        if byte_len % size_of::<PackedReactionRecord>() != 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "sigpack byte length {byte_len} is not divisible by record size {}",
                    size_of::<PackedReactionRecord>()
                ),
            ));
        }

        let map = if byte_len == 0 {
            None
        } else {
            // SAFETY: the file is opened read-only, its length has been
            // validated, and `Mmap` owns the mapping independently of `file`.
            Some(unsafe { MmapOptions::new().map(&file)? })
        };
        Ok(Self { map })
    }

    /// Borrows all mapped records without allocation or deserialization.
    #[must_use]
    pub fn records(&self) -> &[PackedReactionRecord] {
        match &self.map {
            Some(map) => bytemuck::cast_slice(map.as_ref()),
            None => &[],
        }
    }

    /// Number of records in the mapped matrix.
    #[must_use]
    pub fn len(&self) -> usize {
        self.records().len()
    }

    /// Whether the mapped matrix contains no records.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.map.is_none()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::OpenOptions;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temporary_path(test_name: &str) -> std::path::PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!(
            "steric_x_{test_name}_{}_{nonce}.sigpack",
            std::process::id()
        ))
    }

    #[test]
    fn supports_empty_matrices() {
        let path = temporary_path("empty");
        SigPackWriter::export(&[], &path).unwrap();
        let reader = SigPackReader::open(&path).unwrap();
        assert!(reader.is_empty());
        drop(reader);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn rejects_partial_records() {
        let path = temporary_path("partial");
        let mut file = OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(&path)
            .unwrap();
        file.write_all(&[0_u8; 63]).unwrap();
        drop(file);

        let error = SigPackReader::open(&path).err().unwrap();
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn version_two_round_trip_is_zero_copy() {
        let path = temporary_path("v2");
        let mut records = vec![PackedReactionRecordV2::default(); 1_000];
        for (index, record) in records.iter_mut().enumerate() {
            record.reaction.l = index as f32 * 0.25;
            record.buried_volume.max_delta_qvbur_min = index as f32 * 0.5;
        }
        SigPackV2Writer::export(&records, &path).unwrap();
        let reader = SigPackV2Reader::open(&path).unwrap();
        assert_eq!(reader.len(), records.len());
        assert_eq!(reader.records(), records);
        assert_eq!(reader.header().record_count, 1_000);
        drop(reader);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn version_two_reader_rejects_legacy_files() {
        let path = temporary_path("v2_rejects_v1");
        SigPackWriter::export(&[PackedReactionRecord::default()], &path).unwrap();
        let error = SigPackV2Reader::open(&path).err().unwrap();
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        std::fs::remove_file(path).unwrap();
    }
}
