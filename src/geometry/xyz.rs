use glam::Vec3;
use std::error::Error;
use std::fmt::{Display, Formatter};
use std::fs;
use std::io;
use std::path::Path;

/// A single atom with Cartesian position and a Bondi-style van der Waals radius.
#[derive(Clone, Debug, PartialEq)]
pub struct Atom {
    /// Element symbol, normalized only by trimming surrounding whitespace.
    pub element: String,
    /// Cartesian coordinate in ångströms.
    pub position: Vec3,
    /// Van der Waals radius in ångströms.
    pub vdw_radius: f32,
}

impl Atom {
    /// Creates an atom and assigns a conventional van der Waals radius.
    #[must_use]
    pub fn new(element: impl Into<String>, position: Vec3) -> Self {
        let element = element.into();
        Self {
            vdw_radius: van_der_waals_radius(&element),
            element,
            position,
        }
    }
}

/// A molecular structure represented as atoms in source-file order.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Molecule {
    /// Atoms in source-file order.
    pub atoms: Vec<Atom>,
}

impl Molecule {
    /// Loads one molecule from a standard XYZ coordinate file.
    ///
    /// The first line must contain the atom count, the second line is treated
    /// as a comment, and each following atom line must contain
    /// `Element X Y Z`. Coordinates are interpreted in ångströms.
    pub fn from_xyz_file(path: &Path) -> Result<Molecule, Box<dyn Error>> {
        let contents = fs::read_to_string(path)?;
        Ok(parse_xyz(&contents)?)
    }
}

/// Errors produced while loading molecular coordinate files.
#[derive(Debug)]
pub enum GeometryError {
    /// Underlying filesystem error.
    Io(io::Error),
    /// Coordinate data did not satisfy the selected format.
    Format(String),
}

impl Display for GeometryError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(error) => write!(f, "coordinate I/O error: {error}"),
            Self::Format(message) => write!(f, "invalid coordinate data: {message}"),
        }
    }
}

impl Error for GeometryError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Io(error) => Some(error),
            Self::Format(_) => None,
        }
    }
}

impl From<io::Error> for GeometryError {
    fn from(value: io::Error) -> Self {
        Self::Io(value)
    }
}

/// Parses all molecular frames supported by a `.xyz` or `.sdf` file.
pub fn parse_coordinate_file(path: impl AsRef<Path>) -> Result<Vec<Molecule>, GeometryError> {
    let path = path.as_ref();
    let contents = fs::read_to_string(path)?;
    match path
        .extension()
        .and_then(|extension| extension.to_str())
        .map(str::to_ascii_lowercase)
        .as_deref()
    {
        Some("xyz") => Ok(vec![parse_xyz(&contents)?]),
        Some("sdf") | Some("mol") => parse_sdf(&contents),
        extension => Err(GeometryError::Format(format!(
            "unsupported extension {:?} for {}",
            extension,
            path.display()
        ))),
    }
}

/// Parses one standard XYZ frame.
pub fn parse_xyz(input: &str) -> Result<Molecule, GeometryError> {
    let mut lines = input.lines();
    let count_line = lines
        .next()
        .ok_or_else(|| GeometryError::Format("XYZ file is empty".into()))?;
    let atom_count = count_line.trim().parse::<usize>().map_err(|_| {
        GeometryError::Format(format!("invalid XYZ atom count `{}`", count_line.trim()))
    })?;
    lines
        .next()
        .ok_or_else(|| GeometryError::Format("XYZ file has no comment line".into()))?;
    let mut atoms = Vec::with_capacity(atom_count);

    for atom_index in 0..atom_count {
        let line = lines.next().ok_or_else(|| {
            GeometryError::Format(format!(
                "XYZ declares {atom_count} atoms but ends at atom {atom_index}"
            ))
        })?;
        atoms.push(parse_atom_fields(line, atom_index + 3, "XYZ")?);
    }

    Ok(Molecule { atoms })
}

/// Parses one or more V2000 mol blocks from an SDF file.
pub fn parse_sdf(input: &str) -> Result<Vec<Molecule>, GeometryError> {
    let mut molecules = Vec::new();
    for (block_index, block) in input.split("$$$$").enumerate() {
        if block.trim().is_empty() {
            continue;
        }
        // The first record begins at the file start, where an empty first line
        // is a legitimate (blank) molecule title. Every later record begins with
        // the single line terminator that followed the `$$$$` delimiter, which is
        // not part of the record and must be dropped so the fixed four-line
        // header (title, program, comment, counts) stays aligned. Blank titles
        // are common — OpenBabel and RDKit both emit them.
        let block = if block_index == 0 {
            block
        } else {
            block
                .strip_prefix("\r\n")
                .or_else(|| block.strip_prefix('\n'))
                .unwrap_or(block)
        };
        molecules.push(parse_mol_block(block, block_index + 1)?);
    }
    if molecules.is_empty() {
        return Err(GeometryError::Format(
            "SDF contains no molecular records".into(),
        ));
    }
    Ok(molecules)
}

fn parse_mol_block(block: &str, block_number: usize) -> Result<Molecule, GeometryError> {
    // The record separator has already been stripped by `parse_sdf`; keep the
    // header intact here (line 0 is the title, which may be blank).
    let lines: Vec<&str> = block.lines().collect();
    if lines.len() < 4 {
        return Err(GeometryError::Format(format!(
            "SDF record {block_number} has no counts line"
        )));
    }
    let counts = lines[3];
    if counts.contains("V3000") {
        return Err(GeometryError::Format(format!(
            "SDF record {block_number} uses unsupported V3000 syntax"
        )));
    }
    let atom_count = counts
        .get(0..3)
        .unwrap_or(counts)
        .trim()
        .parse::<usize>()
        .or_else(|_| counts.split_whitespace().next().unwrap_or_default().parse())
        .map_err(|_| {
            GeometryError::Format(format!(
                "SDF record {block_number} has an invalid counts line"
            ))
        })?;
    if lines.len() < 4 + atom_count {
        return Err(GeometryError::Format(format!(
            "SDF record {block_number} declares {atom_count} atoms but is truncated"
        )));
    }

    let atoms = lines[4..4 + atom_count]
        .iter()
        .enumerate()
        .map(|(index, line)| parse_atom_fields(line, index + 5, "SDF"))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Molecule { atoms })
}

fn parse_atom_fields(line: &str, line_number: usize, format: &str) -> Result<Atom, GeometryError> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    if fields.len() < 4 {
        return Err(GeometryError::Format(format!(
            "{format} atom line {line_number} has fewer than four fields"
        )));
    }
    let parse_coordinate = |field: &str, axis: &str| {
        field.parse::<f32>().map_err(|_| {
            GeometryError::Format(format!(
                "{format} atom line {line_number} has invalid {axis} coordinate `{field}`"
            ))
        })
    };
    // XYZ uses element,x,y,z while V2000 uses x,y,z,element.
    if fields[0].parse::<f32>().is_ok() {
        let position = Vec3::new(
            parse_coordinate(fields[0], "x")?,
            parse_coordinate(fields[1], "y")?,
            parse_coordinate(fields[2], "z")?,
        );
        Ok(Atom::new(fields[3], position))
    } else {
        let position = Vec3::new(
            parse_coordinate(fields[1], "x")?,
            parse_coordinate(fields[2], "y")?,
            parse_coordinate(fields[3], "z")?,
        );
        Ok(Atom::new(fields[0], position))
    }
}

/// Returns a conventional van der Waals radius in ångströms.
///
/// Unknown elements use a conservative 1.80 Å fallback instead of preventing
/// ingestion of otherwise valid coordinate data.
#[must_use]
pub fn van_der_waals_radius(element: &str) -> f32 {
    match element.trim().to_ascii_uppercase().as_str() {
        "H" => 1.20,
        "B" => 1.92,
        "C" => 1.70,
        "N" => 1.55,
        "O" => 1.52,
        "F" => 1.47,
        "SI" => 2.10,
        "P" => 1.80,
        "S" => 1.80,
        "CL" => 1.75,
        "BR" => 1.85,
        "I" => 1.98,
        _ => 1.80,
    }
}

/// Returns a Cordero (2008) single-bond covalent radius in ångströms.
///
/// Used to decide bonded connectivity from Cartesian geometry when explicit
/// bond records are absent or unreliable. Two atoms are treated as bonded when
/// their separation is within a small tolerance of the summed covalent radii.
/// Unknown elements fall back to 0.77 Å (a carbon-like default) rather than
/// rejecting otherwise valid coordinate data.
#[must_use]
pub fn covalent_radius(element: &str) -> f32 {
    match element.trim().to_ascii_uppercase().as_str() {
        "H" => 0.31,
        "B" => 0.84,
        "C" => 0.76,
        "N" => 0.71,
        "O" => 0.66,
        "F" => 0.57,
        "SI" => 1.11,
        "P" => 1.07,
        "S" => 1.05,
        "CL" => 1.02,
        "BR" => 1.20,
        "I" => 1.39,
        _ => 0.77,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn parses_multi_atom_xyz_string() {
        let sample = "\
4
amine fragment
C  0.000000  0.000000  0.000000
H  1.0e0      0.000000  0.000000
N -1.250000   0.500000  0.000000
O  0.000000  -1.400000  0.250000
";
        let molecule = parse_xyz(sample).unwrap();

        assert_eq!(molecule.atoms.len(), 4);
        assert_eq!(molecule.atoms[0].element, "C");
        assert_eq!(molecule.atoms[0].vdw_radius, 1.70);
        assert_eq!(molecule.atoms[1].position, Vec3::new(1.0, 0.0, 0.0));
        assert_eq!(molecule.atoms[1].vdw_radius, 1.20);
        assert_eq!(molecule.atoms[2].position, Vec3::new(-1.25, 0.5, 0.0));
        assert_eq!(molecule.atoms[2].vdw_radius, 1.55);
        assert_eq!(molecule.atoms[3].position.z, 0.25);
        assert_eq!(molecule.atoms[3].vdw_radius, 1.52);
    }

    #[test]
    fn loads_molecule_from_xyz_file() {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path =
            std::env::temp_dir().join(format!("steric_x_xyz_{}_{nonce}.xyz", std::process::id()));
        fs::write(&path, "3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n").unwrap();

        let molecule = Molecule::from_xyz_file(&path).unwrap();

        assert_eq!(molecule.atoms.len(), 3);
        assert_eq!(molecule.atoms[0].element, "O");
        assert_eq!(molecule.atoms[0].vdw_radius, 1.52);
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn parses_v2000_sdf() {
        let sdf = "water\n  stericx\n\n  3  2  0  0  0  0            999 V2000\n\
                   0.0000    0.0000    0.0000 O   0  0\n\
                   0.9500    0.0000    0.0000 H   0  0\n\
                  -0.2500    0.9200    0.0000 H   0  0\nM  END\n$$$$\n";
        let molecules = parse_sdf(sdf).unwrap();
        assert_eq!(molecules.len(), 1);
        assert_eq!(molecules[0].atoms[0].element, "O");
    }

    #[test]
    fn parses_sdf_with_blank_title_line() {
        // OpenBabel and RDKit routinely write an empty first (title) line.
        let sdf = "\n OpenBabel\n\n  2  1  0  0  0  0            999 V2000\n\
                   0.0000    0.0000    0.0000 P   0  0\n\
                   0.0000    0.0000    1.4300 C   0  0\nM  END\n$$$$\n";
        let molecules = parse_sdf(sdf).unwrap();
        assert_eq!(molecules.len(), 1);
        assert_eq!(molecules[0].atoms.len(), 2);
        assert_eq!(molecules[0].atoms[0].element, "P");
    }

    #[test]
    fn parses_multi_record_sdf_with_blank_titles() {
        // Each record has a blank title; the `$$$$` delimiter sits on its own
        // line, so the newline that follows it must not be mistaken for a title.
        let record = "\n prog\n\n  1  0  0  0  0  0            999 V2000\n\
                      0.0000    0.0000    0.0000 P   0  0\nM  END\n";
        let sdf = format!("{record}$$$$\n{record}$$$$\n");
        let molecules = parse_sdf(&sdf).unwrap();
        assert_eq!(molecules.len(), 2);
        assert!(molecules.iter().all(|m| m.atoms[0].element == "P"));
    }

    #[test]
    fn rejects_truncated_xyz() {
        let error = parse_xyz("2\nname\nC 0 0 0\n").unwrap_err();
        assert!(error.to_string().contains("ends at atom 1"));
    }
}
