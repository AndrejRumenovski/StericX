use bytemuck::{Pod, Zeroable};

/// One cache-line-sized physical-organic reaction observation.
///
/// The layout is exactly sixteen contiguous `f32` values (64 bytes), permitting
/// direct conversion between record slices and native-endian byte slices.
#[repr(C, align(64))]
#[derive(Clone, Copy, Debug, Default, PartialEq, Pod, Zeroable)]
pub struct PackedReactionRecord {
    /// Sterimol length along the primary bond axis, in ångströms.
    pub l: f32,
    /// Minimum perpendicular Sterimol width, in ångströms.
    pub b1: f32,
    /// Maximum perpendicular Sterimol width, in ångströms.
    pub b5: f32,
    /// Donor-atom NBO partial charge.
    pub nbo_charge: f32,
    /// Diagnostic IR stretching frequency, in cm⁻¹.
    pub ir_freq: f32,
    /// Experimental or simulation temperature, in kelvin.
    pub temp_k: f32,
    /// Experimental transition-state energy difference, in kcal/mol.
    pub exp_ddg: f32,
    /// Ensemble metadata preserving the 64-byte binary ABI.
    ///
    /// The slots contain `[L_min, L_max, B1_min, B1_max, B5_min, B5_max,
    /// conformer_count, energy_span_kcal_mol, ensemble_schema_version]`.
    pub reserved: [f32; 9],
}

/// Kraken-style conformer-ensemble buried-volume descriptor block.
///
/// Sixteen contiguous floats retain cache-line alignment and make descriptor
/// matrices suitable for direct SIMD loads.
#[repr(C, align(64))]
#[derive(Clone, Copy, Debug, Default, PartialEq, Pod, Zeroable)]
pub struct PackedBuriedVolumeRecord {
    pub vbur_boltz: f32,
    pub vbur_min: f32,
    pub vbur_max: f32,
    pub vbur_delta: f32,
    pub qvbur_min_boltz: f32,
    pub qvbur_max_boltz: f32,
    pub max_delta_qvbur_boltz: f32,
    pub max_delta_qvbur_min: f32,
    pub max_delta_qvbur_max: f32,
    pub max_delta_qvbur_delta: f32,
    pub max_delta_qvbur_vburminconf: f32,
    pub near_vbur_boltz: f32,
    pub far_vbur_boltz: f32,
    pub conformer_count: f32,
    pub sphere_radius: f32,
    pub grid_density: f32,
}

/// Version-two reaction record containing the stable v1 ABI plus one
/// coordination-aware descriptor cache line.
#[repr(C, align(64))]
#[derive(Clone, Copy, Debug, Default, PartialEq, Pod, Zeroable)]
pub struct PackedReactionRecordV2 {
    pub reaction: PackedReactionRecord,
    pub buried_volume: PackedBuriedVolumeRecord,
}

/// Fixed cache-line header distinguishing v2 files from legacy flat v1 files.
#[repr(C, align(64))]
#[derive(Clone, Copy, Debug, PartialEq, Pod, Zeroable)]
pub struct SigPackHeaderV2 {
    pub magic: [u8; 8],
    pub schema_version: u32,
    pub endian_marker: u32,
    pub record_count: u64,
    pub record_size: u32,
    pub descriptor_count: u32,
    pub reserved: [u8; 32],
}

impl SigPackHeaderV2 {
    pub const MAGIC: [u8; 8] = *b"SIGPKV2\0";
    pub const SCHEMA_VERSION: u32 = 2;
    pub const ENDIAN_MARKER: u32 = 0x0102_0304;
    pub const DESCRIPTOR_COUNT: u32 = 32;

    #[must_use]
    pub fn new(record_count: usize) -> Self {
        Self {
            magic: Self::MAGIC,
            schema_version: Self::SCHEMA_VERSION,
            endian_marker: Self::ENDIAN_MARKER,
            record_count: record_count as u64,
            record_size: size_of::<PackedReactionRecordV2>() as u32,
            descriptor_count: Self::DESCRIPTOR_COUNT,
            reserved: [0; 32],
        }
    }
}

impl PackedReactionRecord {
    /// Current meaning of the nine ensemble metadata slots.
    pub const ENSEMBLE_SCHEMA_VERSION: f32 = 1.0;

    /// Constructs a geometry-only record with remaining fields zeroed.
    #[must_use]
    pub fn from_sterimol(l: f32, b1: f32, b5: f32) -> Self {
        Self {
            l,
            b1,
            b5,
            ..Self::default()
        }
    }

    /// Constructs a Boltzmann-averaged record with conformer envelope metadata.
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn from_ensemble(
        l: f32,
        b1: f32,
        b5: f32,
        l_min: f32,
        l_max: f32,
        b1_min: f32,
        b1_max: f32,
        b5_min: f32,
        b5_max: f32,
        conformer_count: usize,
        energy_span_kcal_mol: f32,
    ) -> Self {
        Self {
            l,
            b1,
            b5,
            reserved: [
                l_min,
                l_max,
                b1_min,
                b1_max,
                b5_min,
                b5_max,
                conformer_count as f32,
                energy_span_kcal_mol,
                Self::ENSEMBLE_SCHEMA_VERSION,
            ],
            ..Self::default()
        }
    }

    /// Minimum conformer Sterimol length.
    #[must_use]
    pub const fn l_min(&self) -> f32 {
        self.reserved[0]
    }

    /// Maximum conformer Sterimol length.
    #[must_use]
    pub const fn l_max(&self) -> f32 {
        self.reserved[1]
    }

    /// Minimum conformer B1.
    #[must_use]
    pub const fn b1_min(&self) -> f32 {
        self.reserved[2]
    }

    /// Maximum conformer B1.
    #[must_use]
    pub const fn b1_max(&self) -> f32 {
        self.reserved[3]
    }

    /// Minimum conformer B5.
    #[must_use]
    pub const fn b5_min(&self) -> f32 {
        self.reserved[4]
    }

    /// Maximum conformer B5.
    #[must_use]
    pub const fn b5_max(&self) -> f32 {
        self.reserved[5]
    }

    /// Number of conformers represented by this record.
    #[must_use]
    pub fn conformer_count(&self) -> usize {
        self.reserved[6].max(0.0).round() as usize
    }

    /// Energy span of the retained ensemble in kcal/mol.
    #[must_use]
    pub const fn ensemble_energy_span(&self) -> f32 {
        self.reserved[7]
    }
}

const _: () = assert!(size_of::<PackedReactionRecord>() == 64);
const _: () = assert!(align_of::<PackedReactionRecord>() == 64);
const _: () = assert!(size_of::<PackedBuriedVolumeRecord>() == 64);
const _: () = assert!(align_of::<PackedBuriedVolumeRecord>() == 64);
const _: () = assert!(size_of::<PackedReactionRecordV2>() == 128);
const _: () = assert!(align_of::<PackedReactionRecordV2>() == 64);
const _: () = assert!(size_of::<SigPackHeaderV2>() == 64);
const _: () = assert!(align_of::<SigPackHeaderV2>() == 64);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn record_is_exactly_one_cache_line_and_pod_castable() {
        let record = PackedReactionRecord::from_sterimol(2.0, 1.0, 3.0);
        assert_eq!(size_of::<PackedReactionRecord>(), 64);
        assert_eq!(align_of::<PackedReactionRecord>(), 64);
        assert_eq!(bytemuck::bytes_of(&record).len(), 64);
    }

    #[test]
    fn ensemble_metadata_round_trips_through_reserved_slots() {
        let record = PackedReactionRecord::from_ensemble(
            2.0, 1.5, 3.0, 1.8, 2.2, 1.4, 1.7, 2.8, 3.3, 12, 4.2,
        );
        assert_eq!(record.conformer_count(), 12);
        assert_eq!(record.l_min(), 1.8);
        assert_eq!(record.b5_max(), 3.3);
        assert_eq!(record.ensemble_energy_span(), 4.2);
        assert_eq!(
            record.reserved[8],
            PackedReactionRecord::ENSEMBLE_SCHEMA_VERSION
        );
    }

    #[test]
    fn version_two_layout_is_header_plus_two_cache_lines_per_record() {
        let header = SigPackHeaderV2::new(1_000);
        assert_eq!(header.magic, SigPackHeaderV2::MAGIC);
        assert_eq!(header.record_size, 128);
        assert_eq!(header.record_count, 1_000);
        assert_eq!(bytemuck::bytes_of(&header).len(), 64);
        assert_eq!(size_of::<PackedReactionRecordV2>(), 128);
    }
}
