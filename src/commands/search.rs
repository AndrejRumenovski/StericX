//! `stericx search`: rank a ligand library by steric similarity to a query,
//! under optional descriptor constraints.
//!
//! Similarity is a Euclidean distance in *standardized* descriptor space: each
//! feature is z-scored against the library's own mean and standard deviation
//! before the distance is taken, so descriptors on different scales (a Sterimol
//! `L` near 8 Å, a `%Vbur` near 30, a `pyr_P` near 0.9) contribute comparably
//! instead of being dominated by whichever happens to carry the largest units.
//! The distance is a *shape-similarity ranking*, not a prediction of reactivity.

use crate::cli::{DescriptorFormat, SterimolAxis};
use crate::descriptors::{DescriptorResult, descriptors_for_file};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::error::Error;
use std::path::{Path, PathBuf};
use steric_x::BuriedVolumeConfig;

/// One library member: the descriptor columns emitted by `stericx descriptors
/// --format csv`, so a CSV produced by that command is a valid library as-is.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub(crate) struct LibraryEntry {
    pub(crate) file: String,
    pub(crate) conformers: usize,
    pub(crate) donor_element: String,
    pub(crate) donor_index: usize,
    pub(crate) substituents: String,
    pub(crate) sterimol_l: f32,
    pub(crate) sterimol_b1: f32,
    pub(crate) sterimol_b5: f32,
    pub(crate) percent_buried_volume: f32,
    pub(crate) buried_volume: f32,
    pub(crate) qvbur_min: f32,
    pub(crate) qvbur_max: f32,
    pub(crate) max_delta_qvbur: f32,
    pub(crate) max_delta_qvbur_min: f32,
    pub(crate) pyr_p: f32,
    pub(crate) pyr_alpha: f32,
}

impl From<&DescriptorResult> for LibraryEntry {
    fn from(result: &DescriptorResult) -> Self {
        Self {
            file: result.file.clone(),
            conformers: result.conformers,
            donor_element: result.donor_element.clone(),
            donor_index: result.donor_index,
            substituents: result.substituents.join(" "),
            sterimol_l: result.sterimol_l,
            sterimol_b1: result.sterimol_b1,
            sterimol_b5: result.sterimol_b5,
            percent_buried_volume: result.percent_buried_volume,
            buried_volume: result.buried_volume,
            qvbur_min: result.qvbur_min,
            qvbur_max: result.qvbur_max,
            max_delta_qvbur: result.max_delta_qvbur,
            max_delta_qvbur_min: result.max_delta_qvbur_min,
            pyr_p: result.pyr_p,
            pyr_alpha: result.pyr_alpha,
        }
    }
}

/// A searchable descriptor: its canonical name, short aliases accepted on the
/// command line, and how to read it off an entry.
struct Feature {
    name: &'static str,
    aliases: &'static [&'static str],
    get: fn(&LibraryEntry) -> f32,
}

const FEATURES: &[Feature] = &[
    Feature {
        name: "sterimol_l",
        aliases: &["l"],
        get: |e| e.sterimol_l,
    },
    Feature {
        name: "sterimol_b1",
        aliases: &["b1"],
        get: |e| e.sterimol_b1,
    },
    Feature {
        name: "sterimol_b5",
        aliases: &["b5"],
        get: |e| e.sterimol_b5,
    },
    Feature {
        name: "percent_buried_volume",
        aliases: &["vbur", "percent_vbur", "%vbur"],
        get: |e| e.percent_buried_volume,
    },
    Feature {
        name: "buried_volume",
        aliases: &["vbur_angstrom3"],
        get: |e| e.buried_volume,
    },
    Feature {
        name: "qvbur_min",
        aliases: &[],
        get: |e| e.qvbur_min,
    },
    Feature {
        name: "qvbur_max",
        aliases: &[],
        get: |e| e.qvbur_max,
    },
    Feature {
        name: "max_delta_qvbur",
        aliases: &["quadrant_asymmetry"],
        get: |e| e.max_delta_qvbur,
    },
    Feature {
        name: "max_delta_qvbur_min",
        aliases: &[],
        get: |e| e.max_delta_qvbur_min,
    },
    Feature {
        name: "pyr_p",
        aliases: &["pyr"],
        get: |e| e.pyr_p,
    },
    Feature {
        name: "pyr_alpha",
        aliases: &[],
        get: |e| e.pyr_alpha,
    },
];

/// The default similarity space: the ligand's shape envelope (`L`/`B1`/`B5`),
/// how much of the metal's coordination sphere it fills (`%Vbur`), how unevenly
/// it fills it (`max_delta_qvbur`), and the donor's pyramidalization.
///
/// `buried_volume` is deliberately absent: it is `percent_buried_volume`
/// rescaled by a constant sphere volume, so including both would silently
/// double-weight the same physical quantity.
const DEFAULT_FEATURES: &[&str] = &[
    "sterimol_l",
    "sterimol_b1",
    "sterimol_b5",
    "percent_buried_volume",
    "max_delta_qvbur",
    "pyr_p",
];

fn resolve_feature(name: &str) -> Result<usize, String> {
    let wanted = name.trim().to_ascii_lowercase();
    FEATURES
        .iter()
        .position(|feature| {
            feature.name == wanted || feature.aliases.iter().any(|alias| *alias == wanted)
        })
        .ok_or_else(|| {
            format!(
                "unknown descriptor `{name}`; available: {}",
                FEATURES
                    .iter()
                    .map(|feature| feature.name)
                    .collect::<Vec<_>>()
                    .join(", ")
            )
        })
}

/// A numeric constraint on one descriptor.
#[derive(Clone, Copy, Debug)]
struct Filter {
    feature: usize,
    bound: Bound,
}

#[derive(Clone, Copy, Debug)]
enum Bound {
    Range(f32, f32),
    Below(f32, bool),
    Above(f32, bool),
}

impl Filter {
    fn accepts(&self, entry: &LibraryEntry) -> bool {
        let value = (FEATURES[self.feature].get)(entry);
        match self.bound {
            Bound::Range(low, high) => value >= low && value <= high,
            Bound::Below(limit, inclusive) => {
                if inclusive {
                    value <= limit
                } else {
                    value < limit
                }
            }
            Bound::Above(limit, inclusive) => {
                if inclusive {
                    value >= limit
                } else {
                    value > limit
                }
            }
        }
    }

    fn describe(&self) -> String {
        let name = FEATURES[self.feature].name;
        match self.bound {
            Bound::Range(low, high) => format!("{name} in [{low}, {high}]"),
            Bound::Below(limit, true) => format!("{name} <= {limit}"),
            Bound::Below(limit, false) => format!("{name} < {limit}"),
            Bound::Above(limit, true) => format!("{name} >= {limit}"),
            Bound::Above(limit, false) => format!("{name} > {limit}"),
        }
    }
}

fn parse_number(text: &str, expression: &str) -> Result<f32, String> {
    text.trim()
        .parse::<f32>()
        .map_err(|_| {
            format!(
                "`{expression}` contains an invalid number `{}`",
                text.trim()
            )
        })
        .and_then(|value| {
            value
                .is_finite()
                .then_some(value)
                .ok_or_else(|| format!("`{expression}` bound must be finite"))
        })
}

/// Parse one `--filter` expression: `name=LOW..HIGH`, or `name` followed by
/// `<`, `<=`, `>`, or `>=` and a value.
fn parse_filter(expression: &str) -> Result<Filter, String> {
    for (token, inclusive, is_below) in [
        ("<=", true, true),
        (">=", true, false),
        ("<", false, true),
        (">", false, false),
    ] {
        if let Some((name, value)) = expression.split_once(token) {
            let feature = resolve_feature(name)?;
            let limit = parse_number(value, expression)?;
            let bound = if is_below {
                Bound::Below(limit, inclusive)
            } else {
                Bound::Above(limit, inclusive)
            };
            return Ok(Filter { feature, bound });
        }
    }

    let (name, range) = expression.split_once('=').ok_or_else(|| {
        format!("`{expression}` is not a filter; use name=LOW..HIGH, name<VALUE, or name>VALUE")
    })?;
    let (low, high) = range
        .split_once("..")
        .ok_or_else(|| format!("`{expression}` must express a range as LOW..HIGH"))?;
    let feature = resolve_feature(name)?;
    let low = parse_number(low, expression)?;
    let high = parse_number(high, expression)?;
    if low > high {
        return Err(format!(
            "`{expression}` has a low bound above its high bound"
        ));
    }
    Ok(Filter {
        feature,
        bound: Bound::Range(low, high),
    })
}

/// Collect coordinate files below a directory, or read a precomputed CSV.
fn load_library(
    library: &Path,
    donor_element: &str,
    sterimol_axis: SterimolAxis,
    config: BuriedVolumeConfig,
) -> Result<Vec<LibraryEntry>, Box<dyn Error>> {
    if library.is_dir() {
        let mut paths = Vec::new();
        collect_coordinate_files(library, &mut paths)?;
        paths.sort();
        if paths.is_empty() {
            return Err(format!(
                "no .xyz/.sdf/.mol geometries found below {}",
                library.display()
            )
            .into());
        }
        // Featurizing a library is embarrassingly parallel and, unlike the
        // `descriptors` command, no committed benchmark measures it single-core.
        let entries = paths
            .par_iter()
            .filter_map(|path| {
                match descriptors_for_file(path, donor_element, None, sterimol_axis, config) {
                    Ok(result) => Some(LibraryEntry::from(&result)),
                    Err(message) => {
                        eprintln!("skipped {}: {message}", path.display());
                        None
                    }
                }
            })
            .collect::<Vec<_>>();
        if entries.is_empty() {
            return Err("no library geometry could be featurized".into());
        }
        return Ok(entries);
    }

    if !library.is_file() {
        return Err(format!("library does not exist: {}", library.display()).into());
    }
    let mut reader = csv::ReaderBuilder::new()
        .trim(csv::Trim::All)
        .from_path(library)?;
    let entries = reader
        .deserialize::<LibraryEntry>()
        .enumerate()
        .map(|(index, row)| {
            row.map_err(|error| {
                format!(
                    "{} row {} could not be parsed: {error}",
                    library.display(),
                    index + 2
                )
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    if entries.is_empty() {
        return Err(format!("library CSV contains no rows: {}", library.display()).into());
    }
    Ok(entries)
}

fn collect_coordinate_files(
    directory: &Path,
    found: &mut Vec<PathBuf>,
) -> Result<(), Box<dyn Error>> {
    for entry in std::fs::read_dir(directory)? {
        let path = entry?.path();
        if path.is_dir() {
            collect_coordinate_files(&path, found)?;
        } else if path
            .extension()
            .and_then(|extension| extension.to_str())
            .is_some_and(|extension| {
                matches!(
                    extension.to_ascii_lowercase().as_str(),
                    "xyz" | "sdf" | "mol"
                )
            })
        {
            found.push(path);
        }
    }
    Ok(())
}

/// Mean and standard deviation of each selected feature across the library.
struct Standardizer {
    features: Vec<usize>,
    means: Vec<f32>,
    deviations: Vec<f32>,
}

impl Standardizer {
    fn fit(entries: &[LibraryEntry], features: &[usize]) -> Self {
        let count = entries.len() as f32;
        let mut means = Vec::with_capacity(features.len());
        let mut deviations = Vec::with_capacity(features.len());
        for &feature in features {
            let read = FEATURES[feature].get;
            let mean = entries.iter().map(read).sum::<f32>() / count;
            let variance = entries
                .iter()
                .map(|entry| (read(entry) - mean).powi(2))
                .sum::<f32>()
                / count;
            means.push(mean);
            deviations.push(variance.sqrt());
        }
        Self {
            features: features.to_vec(),
            means,
            deviations,
        }
    }

    /// Features the library cannot discriminate on (every member identical) are
    /// reported so a ranking is never silently based on fewer axes than asked.
    fn degenerate(&self) -> Vec<&'static str> {
        self.features
            .iter()
            .zip(&self.deviations)
            .filter(|(_, deviation)| **deviation <= f32::EPSILON)
            .map(|(feature, _)| FEATURES[*feature].name)
            .collect()
    }

    fn distance(&self, query: &LibraryEntry, candidate: &LibraryEntry) -> f32 {
        self.features
            .iter()
            .zip(&self.means)
            .zip(&self.deviations)
            .filter(|(_, deviation)| **deviation > f32::EPSILON)
            .map(|((feature, mean), deviation)| {
                let read = FEATURES[*feature].get;
                let query_z = (read(query) - mean) / deviation;
                let candidate_z = (read(candidate) - mean) / deviation;
                (query_z - candidate_z).powi(2)
            })
            .sum::<f32>()
            .sqrt()
    }
}

#[derive(Debug, Serialize)]
struct SearchHit {
    rank: usize,
    distance: f32,
    #[serde(flatten)]
    entry: LibraryEntry,
}

#[derive(Debug, Serialize)]
struct SearchReport {
    query: LibraryEntry,
    similarity_features: Vec<String>,
    filters: Vec<String>,
    library_size: usize,
    candidates_after_filters: usize,
    hits: Vec<SearchHit>,
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn search_command(
    ligand: &Path,
    library: &Path,
    top: usize,
    feature_names: Option<&str>,
    filter_expressions: &[String],
    less_bulky: bool,
    more_bulky: bool,
    donor_element: &str,
    sterimol_axis: SterimolAxis,
    format: DescriptorFormat,
    config: BuriedVolumeConfig,
) -> Result<(), Box<dyn Error>> {
    if less_bulky && more_bulky {
        return Err("--less-bulky and --more-bulky are mutually exclusive".into());
    }
    if top == 0 {
        return Err("--top must be at least 1".into());
    }

    let features = match feature_names {
        Some(list) => list
            .split(',')
            .map(str::trim)
            .filter(|name| !name.is_empty())
            .map(resolve_feature)
            .collect::<Result<Vec<_>, _>>()?,
        None => DEFAULT_FEATURES
            .iter()
            .map(|name| resolve_feature(name))
            .collect::<Result<Vec<_>, _>>()?,
    };
    if features.is_empty() {
        return Err("--features selected no descriptors".into());
    }

    let mut filters = filter_expressions
        .iter()
        .map(|expression| parse_filter(expression))
        .collect::<Result<Vec<_>, _>>()?;

    let query_result = descriptors_for_file(ligand, donor_element, None, sterimol_axis, config)
        .map_err(|message| format!("query ligand: {message}"))?;
    let query = LibraryEntry::from(&query_result);

    // "Find me a less bulky alternative" is a constraint relative to the query,
    // so it can only be resolved once the query has been featurized.
    let vbur = resolve_feature("percent_buried_volume")?;
    if less_bulky {
        filters.push(Filter {
            feature: vbur,
            bound: Bound::Below(query.percent_buried_volume, false),
        });
    }
    if more_bulky {
        filters.push(Filter {
            feature: vbur,
            bound: Bound::Above(query.percent_buried_volume, false),
        });
    }

    let entries = load_library(library, donor_element, sterimol_axis, config)?;
    let library_size = entries.len();
    if library_size < 2 {
        return Err(
            "a library needs at least two members to standardize descriptors for ranking".into(),
        );
    }

    let standardizer = Standardizer::fit(&entries, &features);
    let degenerate = standardizer.degenerate();
    if degenerate.len() == features.len() {
        return Err(
            "every selected descriptor is constant across the library, so nothing can be ranked"
                .into(),
        );
    }
    for name in &degenerate {
        eprintln!("note: `{name}` is constant across the library and cannot affect the ranking");
    }

    let candidates = entries
        .iter()
        .filter(|entry| entry.file != query.file)
        .filter(|entry| filters.iter().all(|filter| filter.accepts(entry)))
        .collect::<Vec<_>>();
    let candidates_after_filters = candidates.len();

    let mut ranked = candidates
        .into_iter()
        .map(|entry| (standardizer.distance(&query, entry), entry))
        .collect::<Vec<_>>();
    // Ties break on file name so repeated runs emit an identical ordering.
    ranked.sort_by(|(left, left_entry), (right, right_entry)| {
        left.total_cmp(right)
            .then_with(|| left_entry.file.cmp(&right_entry.file))
    });
    ranked.truncate(top);

    let report = SearchReport {
        query: query.clone(),
        similarity_features: features
            .iter()
            .map(|feature| FEATURES[*feature].name.to_owned())
            .collect(),
        filters: filters.iter().map(Filter::describe).collect(),
        library_size,
        candidates_after_filters,
        hits: ranked
            .into_iter()
            .enumerate()
            .map(|(index, (distance, entry))| SearchHit {
                rank: index + 1,
                distance,
                entry: entry.clone(),
            })
            .collect(),
    };

    match format {
        DescriptorFormat::Text => print_search_text(&report),
        DescriptorFormat::Json => println!("{}", serde_json::to_string_pretty(&report)?),
        DescriptorFormat::Csv => print_search_csv(&report),
    }
    Ok(())
}

fn print_search_text(report: &SearchReport) {
    println!("query          {}", report.query.file);
    println!(
        "               L {:.2}   B1 {:.2}   B5 {:.2}   Vbur {:.1}%   max_delta_qvbur {:.2}   pyr_P {:.3}",
        report.query.sterimol_l,
        report.query.sterimol_b1,
        report.query.sterimol_b5,
        report.query.percent_buried_volume,
        report.query.max_delta_qvbur,
        report.query.pyr_p
    );
    println!("similarity     {}", report.similarity_features.join(", "));
    if report.filters.is_empty() {
        println!("filters        none");
    } else {
        println!("filters        {}", report.filters.join("; "));
    }
    println!(
        "library        {} members, {} passed the filters",
        report.library_size, report.candidates_after_filters
    );
    if report.hits.is_empty() {
        println!("\nno candidate satisfied every constraint");
        return;
    }
    println!(
        "\n{:>4}  {:>8}  {:>6}  {:>6}  {:>6}  {:>7}  {:>7}  {:>6}  file",
        "rank", "distance", "L", "B1", "B5", "%Vbur", "dqvbur", "pyr_P"
    );
    for hit in &report.hits {
        println!(
            "{:>4}  {:>8.3}  {:>6.2}  {:>6.2}  {:>6.2}  {:>7.1}  {:>7.2}  {:>6.3}  {}",
            hit.rank,
            hit.distance,
            hit.entry.sterimol_l,
            hit.entry.sterimol_b1,
            hit.entry.sterimol_b5,
            hit.entry.percent_buried_volume,
            hit.entry.max_delta_qvbur,
            hit.entry.pyr_p,
            hit.entry.file
        );
    }
}

fn print_search_csv(report: &SearchReport) {
    println!(
        "rank,distance,file,donor_element,sterimol_l,sterimol_b1,sterimol_b5,\
         percent_buried_volume,buried_volume,qvbur_min,qvbur_max,max_delta_qvbur,\
         max_delta_qvbur_min,pyr_p,pyr_alpha"
    );
    for hit in &report.hits {
        let entry = &hit.entry;
        println!(
            "{},{:.6},{},{},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4}",
            hit.rank,
            hit.distance,
            crate::descriptors::csv_field(&entry.file),
            entry.donor_element,
            entry.sterimol_l,
            entry.sterimol_b1,
            entry.sterimol_b5,
            entry.percent_buried_volume,
            entry.buried_volume,
            entry.qvbur_min,
            entry.qvbur_max,
            entry.max_delta_qvbur,
            entry.max_delta_qvbur_min,
            entry.pyr_p,
            entry.pyr_alpha
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn entry(file: &str, l: f32, b1: f32, b5: f32, vbur: f32) -> LibraryEntry {
        LibraryEntry {
            file: file.into(),
            conformers: 1,
            donor_element: "P".into(),
            donor_index: 0,
            substituents: "C C C".into(),
            sterimol_l: l,
            sterimol_b1: b1,
            sterimol_b5: b5,
            percent_buried_volume: vbur,
            buried_volume: vbur * 1.7959,
            qvbur_min: 9.0,
            qvbur_max: 13.0,
            max_delta_qvbur: 4.0,
            max_delta_qvbur_min: 4.0,
            pyr_p: 0.93,
            pyr_alpha: 17.0,
        }
    }

    #[test]
    fn parses_every_supported_filter_form() {
        assert!(matches!(
            parse_filter("percent_buried_volume=30..35").unwrap().bound,
            Bound::Range(low, high) if (low - 30.0).abs() < 1e-6 && (high - 35.0).abs() < 1e-6
        ));
        assert!(matches!(
            parse_filter("b5<7").unwrap().bound,
            Bound::Below(limit, false) if (limit - 7.0).abs() < 1e-6
        ));
        assert!(matches!(
            parse_filter("b5<=7").unwrap().bound,
            Bound::Below(_, true)
        ));
        assert!(matches!(
            parse_filter("l>8").unwrap().bound,
            Bound::Above(_, false)
        ));
        assert!(matches!(
            parse_filter("l>=8").unwrap().bound,
            Bound::Above(_, true)
        ));
    }

    #[test]
    fn resolves_short_aliases_to_canonical_descriptors() {
        assert_eq!(
            resolve_feature("vbur").unwrap(),
            resolve_feature("percent_buried_volume").unwrap()
        );
        assert_eq!(
            resolve_feature("B5").unwrap(),
            resolve_feature("sterimol_b5").unwrap()
        );
        assert!(resolve_feature("not_a_descriptor").is_err());
    }

    #[test]
    fn rejects_malformed_filters() {
        assert!(parse_filter("percent_buried_volume").is_err());
        assert!(parse_filter("percent_buried_volume=30").is_err());
        assert!(parse_filter("percent_buried_volume=35..30").is_err());
        assert!(parse_filter("b5<banana").is_err());
    }

    #[test]
    fn filters_select_the_intended_range() {
        let candidate = entry("a.xyz", 8.0, 1.7, 7.0, 32.0);
        assert!(parse_filter("vbur=30..35").unwrap().accepts(&candidate));
        assert!(!parse_filter("vbur=10..20").unwrap().accepts(&candidate));
        assert!(parse_filter("b5<7.5").unwrap().accepts(&candidate));
        assert!(!parse_filter("b5<7.0").unwrap().accepts(&candidate));
    }

    #[test]
    fn standardized_distance_is_scale_invariant_and_self_zero() {
        let entries = vec![
            entry("a.xyz", 8.0, 1.7, 7.0, 30.0),
            entry("b.xyz", 9.0, 2.0, 8.0, 35.0),
            entry("c.xyz", 7.0, 1.5, 6.0, 25.0),
        ];
        let features = vec![
            resolve_feature("sterimol_l").unwrap(),
            resolve_feature("percent_buried_volume").unwrap(),
        ];
        let standardizer = Standardizer::fit(&entries, &features);
        // A ligand is always its own nearest neighbour.
        assert!(standardizer.distance(&entries[0], &entries[0]).abs() < 1e-6);
        // `b` and `c` sit one step either side of `a`, so both are equidistant
        // even though %Vbur spans five times the range that L does.
        let to_b = standardizer.distance(&entries[0], &entries[1]);
        let to_c = standardizer.distance(&entries[0], &entries[2]);
        assert!((to_b - to_c).abs() < 1e-5, "{to_b} vs {to_c}");
    }

    #[test]
    fn constant_descriptors_are_reported_as_degenerate() {
        let entries = vec![
            entry("a.xyz", 8.0, 1.7, 7.0, 30.0),
            entry("b.xyz", 9.0, 1.7, 8.0, 35.0),
        ];
        let features = vec![
            resolve_feature("sterimol_b1").unwrap(),
            resolve_feature("sterimol_l").unwrap(),
        ];
        let standardizer = Standardizer::fit(&entries, &features);
        assert_eq!(standardizer.degenerate(), vec!["sterimol_b1"]);
    }
}
