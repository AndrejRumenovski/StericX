//! Cache-aligned reaction records and zero-copy `.sigpack` storage.

mod pack;
mod schema;

pub use pack::{SigPackReader, SigPackV2Reader, SigPackV2Writer, SigPackWriter};
pub use schema::{
    PackedBuriedVolumeRecord, PackedReactionRecord, PackedReactionRecordV2, SigPackHeaderV2,
};
