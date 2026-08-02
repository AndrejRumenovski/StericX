use super::expand_features;
use crate::storage::PackedReactionRecord;
use rayon::prelude::*;

/// Eight-coefficient multivariate linear regression predictor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RegressXPredictor {
    weights: [f32; 8],
}

impl RegressXPredictor {
    /// Creates a predictor. `weights[0]` is the intercept coefficient because
    /// `expand_features` starts with a constant `1.0`.
    #[must_use]
    pub const fn new(weights: [f32; 8]) -> Self {
        Self { weights }
    }

    /// Returns the model coefficients in feature order.
    #[must_use]
    pub const fn weights(&self) -> &[f32; 8] {
        &self.weights
    }

    /// Predicts ΔΔG‡ for one packed reaction.
    #[must_use]
    #[inline]
    pub fn predict(&self, record: &PackedReactionRecord) -> f32 {
        let features = expand_features(record);
        dot_product(&features, &self.weights)
    }

    /// Predicts ΔΔG‡ for a record slice in parallel while preserving row order.
    ///
    /// Rayon distributes independent rows across its work-stealing thread pool.
    /// On x86/x86-64 with AVX2, each row is evaluated with one eight-lane SIMD
    /// multiply followed by a horizontal reduction. Other targets use an
    /// explicitly unrolled four-lane chunk kernel suitable for auto-vectorizing.
    #[must_use]
    pub fn predict_batch(&self, records: &[PackedReactionRecord]) -> Vec<f32> {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        if std::is_x86_feature_detected!("avx2") {
            return records
                .par_iter()
                .map(|record| {
                    let features = expand_features(record);
                    // SAFETY: AVX2 support is checked once before the parallel
                    // traversal; both arrays contain exactly eight f32 lanes.
                    unsafe { dot_avx2(&features, &self.weights) }
                })
                .collect();
        }

        records
            .par_iter()
            .map(|record| {
                let features = expand_features(record);
                dot_chunked(&features, &self.weights)
            })
            .collect()
    }
}

#[inline]
fn dot_product(features: &[f32; 8], weights: &[f32; 8]) -> f32 {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    if std::is_x86_feature_detected!("avx2") {
        // SAFETY: guarded by runtime AVX2 feature detection.
        return unsafe { dot_avx2(features, weights) };
    }
    dot_chunked(features, weights)
}

#[inline(always)]
fn dot_chunked(features: &[f32; 8], weights: &[f32; 8]) -> f32 {
    let low = features[0] * weights[0]
        + features[1] * weights[1]
        + features[2] * weights[2]
        + features[3] * weights[3];
    let high = features[4] * weights[4]
        + features[5] * weights[5]
        + features[6] * weights[6]
        + features[7] * weights[7];
    low + high
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn dot_avx2(features: &[f32; 8], weights: &[f32; 8]) -> f32 {
    use core::arch::x86_64::{
        _mm_add_ps, _mm_cvtss_f32, _mm_hadd_ps, _mm256_castps256_ps128, _mm256_extractf128_ps,
        _mm256_loadu_ps, _mm256_mul_ps,
    };

    // SAFETY: each unaligned load reads exactly eight initialized f32 values.
    let products = unsafe {
        _mm256_mul_ps(
            _mm256_loadu_ps(features.as_ptr()),
            _mm256_loadu_ps(weights.as_ptr()),
        )
    };
    let low = _mm256_castps256_ps128(products);
    let high = _mm256_extractf128_ps::<1>(products);
    let pairwise = _mm_hadd_ps(_mm_add_ps(low, high), _mm_add_ps(low, high));
    _mm_cvtss_f32(_mm_hadd_ps(pairwise, pairwise))
}

#[cfg(target_arch = "x86")]
#[target_feature(enable = "avx2")]
unsafe fn dot_avx2(features: &[f32; 8], weights: &[f32; 8]) -> f32 {
    use core::arch::x86::{
        _mm_add_ps, _mm_cvtss_f32, _mm_hadd_ps, _mm256_castps256_ps128, _mm256_extractf128_ps,
        _mm256_loadu_ps, _mm256_mul_ps,
    };

    // SAFETY: each unaligned load reads exactly eight initialized f32 values.
    let products = unsafe {
        _mm256_mul_ps(
            _mm256_loadu_ps(features.as_ptr()),
            _mm256_loadu_ps(weights.as_ptr()),
        )
    };
    let low = _mm256_castps256_ps128(products);
    let high = _mm256_extractf128_ps::<1>(products);
    let pairwise = _mm_hadd_ps(_mm_add_ps(low, high), _mm_add_ps(low, high));
    _mm_cvtss_f32(_mm_hadd_ps(pairwise, pairwise))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Instant;

    #[test]
    fn prediction_matches_known_dot_product() {
        let predictor = RegressXPredictor::new([0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.01]);
        let record = PackedReactionRecord {
            l: 2.0,
            b1: 3.0,
            b5: 4.0,
            nbo_charge: -0.5,
            ir_freq: 1_600.0,
            ..PackedReactionRecord::default()
        };
        let expected = expand_features(&record)
            .iter()
            .zip(predictor.weights())
            .map(|(feature, weight)| feature * weight)
            .sum::<f32>();

        assert!((predictor.predict(&record) - expected).abs() < 1.0e-5);
    }

    #[test]
    fn parallel_batch_preserves_order() {
        let predictor = RegressXPredictor::new([0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let records = [
            PackedReactionRecord::from_sterimol(1.0, 0.0, 0.0),
            PackedReactionRecord::from_sterimol(2.0, 0.0, 0.0),
        ];
        assert_eq!(predictor.predict_batch(&records), vec![1.0, 2.0]);
    }

    #[test]
    fn processes_one_hundred_thousand_records_at_high_throughput() {
        const RECORD_COUNT: usize = 100_000;
        let records: Vec<_> = (0..RECORD_COUNT)
            .map(|index| PackedReactionRecord {
                l: index as f32 * 0.001,
                b1: 1.5,
                b5: 3.5,
                nbo_charge: -0.4,
                ir_freq: 1_650.0,
                ..PackedReactionRecord::default()
            })
            .collect();
        let predictor = RegressXPredictor::new([0.1, 0.2, -0.1, 0.3, 0.5, -0.2, 0.1, 0.001]);

        let started = Instant::now();
        let predictions = predictor.predict_batch(&records);
        let elapsed = started.elapsed();
        let records_per_second = RECORD_COUNT as f64 / elapsed.as_secs_f64().max(f64::EPSILON);

        assert_eq!(predictions.len(), RECORD_COUNT);
        assert!(predictions.iter().all(|prediction| prediction.is_finite()));
        assert!(
            records_per_second > 5_000.0,
            "inference throughput was only {records_per_second:.0} records/s"
        );
    }

    /// The unsafe AVX2 kernel must return the same value as the portable scalar
    /// fallback: both evaluate the identical eight-lane inner product, so any
    /// divergence would mean the SIMD path is silently producing different numbers.
    /// The two sum their terms in different orders, so f32 non-associativity makes
    /// the honest check "agree to tolerance," not "bit-for-bit identical." Runs only
    /// where AVX2 is actually available (calling the kernel otherwise is undefined).
    #[test]
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    fn avx2_matches_scalar_fallback() {
        if !std::is_x86_feature_detected!("avx2") {
            return; // No AVX2 on this host; the scalar fallback is the only path used.
        }
        // Deterministic xorshift so a failure is always reproducible.
        let mut seed: u64 = 0x2545_F491_4F6C_DD1D;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 40) as f32 / (1u32 << 24) as f32 // uniform in [0, 1)
        };
        for _ in 0..5_000 {
            // Spanning realistic descriptor and weight magnitudes, both signs.
            let features: [f32; 8] = std::array::from_fn(|_| (next() - 0.5) * 200.0);
            let weights: [f32; 8] = std::array::from_fn(|_| (next() - 0.5) * 4.0);
            // SAFETY: AVX2 availability is checked at the top of this test.
            let simd = unsafe { dot_avx2(&features, &weights) };
            let scalar = dot_chunked(&features, &weights);
            // Tolerance scaled to the magnitude of the summed terms, so it survives
            // cancellation (a near-zero result from large opposite-sign products).
            let magnitude: f32 = features
                .iter()
                .zip(&weights)
                .map(|(feature, weight)| (feature * weight).abs())
                .sum();
            let tolerance = 1.0e-3 * magnitude + 1.0e-4;
            assert!(
                (simd - scalar).abs() <= tolerance,
                "AVX2 {simd} vs scalar {scalar} diverged beyond tolerance {tolerance}"
            );
        }
    }
}
