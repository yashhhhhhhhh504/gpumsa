use anyhow::{bail, Result};

use crate::fasta::Sequence;
use crate::gpu::score_composition_matrix;

const AA_BUCKETS: usize = 21;

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum BackendKind {
    Auto,
    Cpu,
    Gpu,
}

#[derive(Debug)]
pub struct SimilarityResult {
    pub matrix: Vec<f32>,
    pub backend: &'static str,
}

pub fn compute_similarity_matrix(sequences: &[Sequence], backend: BackendKind) -> Result<SimilarityResult> {
    if sequences.is_empty() {
        bail!("cannot score an empty sequence collection");
    }

    let features = composition_features(sequences);
    match backend {
        BackendKind::Cpu => Ok(SimilarityResult {
            matrix: cpu_similarity(&features, sequences.len(), AA_BUCKETS),
            backend: "cpu",
        }),
        BackendKind::Gpu => Ok(SimilarityResult {
            matrix: score_composition_matrix(sequences.len(), AA_BUCKETS, &features)?,
            backend: "gpu",
        }),
        BackendKind::Auto => match score_composition_matrix(sequences.len(), AA_BUCKETS, &features) {
            Ok(matrix) => Ok(SimilarityResult {
                matrix,
                backend: "gpu",
            }),
            Err(_) => Ok(SimilarityResult {
                matrix: cpu_similarity(&features, sequences.len(), AA_BUCKETS),
                backend: "cpu-fallback",
            }),
        },
    }
}

fn composition_features(sequences: &[Sequence]) -> Vec<f32> {
    let mut features = vec![0.0f32; sequences.len() * AA_BUCKETS];

    for (sequence_index, sequence) in sequences.iter().enumerate() {
        let start = sequence_index * AA_BUCKETS;
        for &residue in &sequence.residues {
            features[start + bucket_for_residue(residue)] += 1.0;
        }

        let norm = features[start..start + AA_BUCKETS]
            .iter()
            .map(|value| value * value)
            .sum::<f32>()
            .sqrt();

        if norm == 0.0 {
            features[start + AA_BUCKETS - 1] = 1.0;
        } else {
            for value in &mut features[start..start + AA_BUCKETS] {
                *value /= norm;
            }
        }
    }

    features
}

fn cpu_similarity(features: &[f32], count: usize, feature_dim: usize) -> Vec<f32> {
    let mut matrix = vec![0.0f32; count * count];

    for i in 0..count {
        for j in 0..count {
            let mut score = 0.0f32;
            let left = i * feature_dim;
            let right = j * feature_dim;
            for offset in 0..feature_dim {
                score += features[left + offset] * features[right + offset];
            }
            matrix[i * count + j] = score;
        }
    }

    matrix
}

fn bucket_for_residue(residue: u8) -> usize {
    match residue {
        b'A' => 0,
        b'C' => 1,
        b'D' => 2,
        b'E' => 3,
        b'F' => 4,
        b'G' => 5,
        b'H' => 6,
        b'I' => 7,
        b'K' => 8,
        b'L' => 9,
        b'M' => 10,
        b'N' => 11,
        b'P' => 12,
        b'Q' => 13,
        b'R' => 14,
        b'S' => 15,
        b'T' => 16,
        b'V' => 17,
        b'W' => 18,
        b'Y' => 19,
        _ => 20,
    }
}

#[cfg(test)]
mod tests {
    use super::{compute_similarity_matrix, BackendKind};
    use crate::fasta::Sequence;

    #[test]
    fn cpu_similarity_is_symmetric() {
        let sequences = vec![
            Sequence {
                id: "a".into(),
                residues: b"AAAA".to_vec(),
            },
            Sequence {
                id: "b".into(),
                residues: b"AAAY".to_vec(),
            },
            Sequence {
                id: "c".into(),
                residues: b"YYYY".to_vec(),
            },
        ];

        let result = compute_similarity_matrix(&sequences, BackendKind::Cpu).unwrap();
        let n = sequences.len();
        for i in 0..n {
            for j in 0..n {
                assert!((result.matrix[i * n + j] - result.matrix[j * n + i]).abs() < 1e-6);
            }
        }
    }
}
