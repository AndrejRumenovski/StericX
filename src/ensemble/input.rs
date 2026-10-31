use super::EnergyMetadata;
use crate::geometry::{MolecularFrame, parse_sdf_with_topology, parse_xyz};
use serde::Serialize;
use std::path::Path;

/// Raw SDF properties are retained in source order, including duplicate keys.
/// Their values are not automatically assigned an energy interpretation.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct SourceProperty {
    pub name: String,
    pub value: String,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Conformer {
    pub index: usize,
    pub identifier: String,
    pub source: String,
    pub source_format: String,
    pub source_geometry: String,
    pub properties: Vec<SourceProperty>,
    pub frame: Result<MolecularFrame, String>,
    pub energy: Option<EnergyMetadata>,
}

/// Ordered supplied structures. File loading does not assert chemical identity
/// across records; callers must establish atom mapping before comparing sectors.
#[derive(Clone, Debug, PartialEq)]
pub struct Ensemble {
    pub conformers: Vec<Conformer>,
}

impl Ensemble {
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self, String> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
        let source = path.display().to_string();
        match path
            .extension()
            .and_then(|x| x.to_str())
            .map(str::to_ascii_lowercase)
            .as_deref()
        {
            Some("sdf" | "mol") => Self::from_sdf_str(&text, &source),
            Some("xyz") => Self::from_xyz_str(&text, &source),
            _ => Err("supported coordinate formats are XYZ and V2000 SDF/MOL".into()),
        }
    }

    pub fn from_sdf(path: impl AsRef<Path>) -> Result<Self, String> {
        let path = path.as_ref();
        Self::from_sdf_str(
            &std::fs::read_to_string(path).map_err(|e| e.to_string())?,
            &path.display().to_string(),
        )
    }

    pub fn from_sdf_str(text: &str, source: &str) -> Result<Self, String> {
        let mut records = Vec::new();
        let mut start = 0;
        let mut end = 0;
        for line in text.split_inclusive('\n') {
            end += line.len();
            if line.trim() == "$$$$" {
                records.push(&text[start..end]);
                start = end;
            }
        }
        if !text[start..].trim().is_empty() {
            records.push(&text[start..]);
        }
        if records.is_empty() {
            return Err("SDF contains no records".into());
        }
        Ok(Self {
            conformers: records
                .into_iter()
                .enumerate()
                .map(|(index, record)| {
                    let frame = parse_sdf_with_topology(record)
                        .map_err(|e| e.to_string())
                        .and_then(|mut frames| {
                            if frames.len() != 1 {
                                return Err("record does not contain exactly one geometry".into());
                            }
                            Ok(frames.remove(0))
                        });
                    Conformer {
                        index,
                        identifier: record.lines().next().unwrap_or("").to_owned(),
                        source: source.into(),
                        source_format: "sdf_v2000".into(),
                        source_geometry: record.into(),
                        properties: properties(record),
                        frame,
                        energy: None,
                    }
                })
                .collect(),
        })
    }

    /// Read concatenated XYZ frames. An invalid atom count makes subsequent
    /// boundaries unknowable: retain the entire remaining source as one error
    /// record, rather than guessing new frame boundaries or dropping text.
    pub fn from_xyz_str(text: &str, source: &str) -> Result<Self, String> {
        let lines: Vec<&str> = text.split_inclusive('\n').collect();
        let mut conformers = Vec::new();
        let mut offset = 0;
        while offset < lines.len() {
            if lines[offset..].iter().all(|l| l.trim().is_empty()) {
                break;
            }
            let count = lines[offset].trim().parse::<usize>();
            let end = count
                .as_ref()
                .ok()
                .and_then(|n| n.checked_add(2))
                .and_then(|n| offset.checked_add(n))
                .unwrap_or(lines.len())
                .min(lines.len());
            let record = lines[offset..end].concat();
            let parsed = parse_xyz(&record)
                .map(|molecule| MolecularFrame {
                    molecule,
                    bonds: None,
                })
                .map_err(|e| e.to_string());
            conformers.push(Conformer {
                index: conformers.len(),
                identifier: lines
                    .get(offset + 1)
                    .unwrap_or(&"")
                    .trim_end_matches(['\r', '\n'])
                    .to_string(),
                source: source.into(),
                source_format: "xyz".into(),
                source_geometry: record,
                properties: Vec::new(),
                frame: parsed,
                energy: None,
            });
            offset = end;
        }
        if conformers.is_empty() {
            return Err("XYZ contains no frames".into());
        }
        Ok(Self { conformers })
    }

    /// Attach declared energy metadata to an original positional record.
    pub fn set_energy(&mut self, index: usize, energy: EnergyMetadata) -> Result<(), String> {
        energy.validate()?;
        let conformer = self
            .conformers
            .get_mut(index)
            .ok_or("conformer index out of bounds")?;
        conformer.energy = Some(energy);
        Ok(())
    }
}

fn properties(record: &str) -> Vec<SourceProperty> {
    let mut result = Vec::new();
    let mut lines = record.lines().peekable();
    // Only property blocks after M END are metadata.
    if !lines.by_ref().any(|line| line.trim() == "M  END") {
        return result;
    }
    while let Some(line) = lines.next() {
        if !line.starts_with('>') {
            continue;
        }
        let Some(name) = line
            .split_once('<')
            .and_then(|(_, rest)| rest.split_once('>'))
            .map(|(name, _)| name)
        else {
            continue;
        };
        let mut value = Vec::new();
        while let Some(line) = lines.peek() {
            if line.trim().is_empty() || line.trim() == "$$$$" || line.starts_with('>') {
                break;
            }
            value.push(lines.next().unwrap());
        }
        result.push(SourceProperty {
            name: name.into(),
            value: value.join("\n"),
        });
    }
    result
}
