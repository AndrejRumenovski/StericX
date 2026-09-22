//! Molecular-coordinate parsing and three-dimensional steric descriptors.

mod buried_volume;
mod pyramidalization;
mod sterimol;
mod xyz;

pub use buried_volume::{
    BuriedVolumeCalculator, BuriedVolumeConfig, BuriedVolumeEnsembleParams, BuriedVolumeError,
    BuriedVolumeParams, bonded_neighbors, coordination_center, coordination_center_with_neighbors,
};
pub use pyramidalization::{PyramidalizationCalculator, PyramidalizationParams};
pub use sterimol::{SterimolCalculator, SterimolParams};
pub use xyz::{
    Atom, GeometryError, MolecularFrame, Molecule, covalent_radius, parse_coordinate_file,
    parse_coordinate_file_with_topology, parse_sdf, parse_sdf_with_topology, parse_xyz,
};

/// Invalid, degenerate, or unrepresentable input to a Cartesian descriptor.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DescriptorError(pub(crate) String);

impl std::fmt::Display for DescriptorError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(&self.0)
    }
}

impl std::error::Error for DescriptorError {}
