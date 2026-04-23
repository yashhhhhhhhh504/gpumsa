use anyhow::{bail, Result};

use crate::fasta::Sequence;
use crate::guide::GuideTree;

#[derive(Copy, Clone, Debug)]
pub struct ScoreScheme {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_penalty: i32,
}

#[derive(Clone, Debug)]
pub struct Alignment {
    pub ids: Vec<String>,
    pub rows: Vec<Vec<u8>>,
}

impl Alignment {
    pub fn from_sequence(sequence: &Sequence) -> Self {
        Self {
            ids: vec![sequence.id.clone()],
            rows: vec![sequence.residues.clone()],
        }
    }

    pub fn column_count(&self) -> usize {
        self.rows.first().map_or(0, Vec::len)
    }

    pub fn without_id(&self, id: &str) -> Option<Self> {
        let keep_indices = self
            .ids
            .iter()
            .enumerate()
            .filter_map(|(index, current)| (current != id).then_some(index))
            .collect::<Vec<_>>();

        if keep_indices.len() == self.ids.len() || keep_indices.is_empty() {
            return None;
        }

        Some(Self {
            ids: keep_indices.iter().map(|&index| self.ids[index].clone()).collect(),
            rows: keep_indices.iter().map(|&index| self.rows[index].clone()).collect(),
        })
    }

    fn profile_tokens(&self) -> Vec<u8> {
        let mut profile = Vec::with_capacity(self.column_count());
        for column in 0..self.column_count() {
            let mut counts = [0usize; 256];
            for row in &self.rows {
                counts[row[column] as usize] += 1;
            }
            let (token, _) = counts
                .iter()
                .enumerate()
                .max_by_key(|(_, count)| **count)
                .expect("non-empty histogram");
            profile.push(token as u8);
        }
        profile
    }

    fn aligned_lengths_match(&self) -> bool {
        let Some(expected) = self.rows.first().map(Vec::len) else {
            return true;
        };
        self.rows.iter().all(|row| row.len() == expected)
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum TraceOp {
    Align,
    GapInLeft,
    GapInRight,
}

pub fn progressive_align(tree: &GuideTree, sequences: &[Sequence], scheme: ScoreScheme) -> Result<Alignment> {
    let alignment = match tree {
        GuideTree::Leaf(index) => Alignment::from_sequence(&sequences[*index]),
        GuideTree::Merge { left, right } => {
            let left_alignment = progressive_align(left, sequences, scheme)?;
            let right_alignment = progressive_align(right, sequences, scheme)?;
            merge_alignments(&left_alignment, &right_alignment, scheme)?
        }
    };

    if !alignment.aligned_lengths_match() {
        bail!("internal alignment invariant failed: rows had inconsistent lengths");
    }

    Ok(alignment)
}

pub fn polish_alignment(
    mut alignment: Alignment,
    sequences: &[Sequence],
    passes: usize,
    scheme: ScoreScheme,
) -> Result<Alignment> {
    if passes == 0 || sequences.len() < 2 {
        return Ok(alignment);
    }

    for _ in 0..passes {
        for sequence in sequences {
            let Some(profile) = alignment.without_id(&sequence.id) else {
                continue;
            };
            alignment = merge_alignments(&profile, &Alignment::from_sequence(sequence), scheme)?;
        }
    }

    Ok(alignment)
}

fn merge_alignments(left: &Alignment, right: &Alignment, scheme: ScoreScheme) -> Result<Alignment> {
    if left.rows.is_empty() || right.rows.is_empty() {
        bail!("cannot merge empty alignments");
    }
    if !left.aligned_lengths_match() || !right.aligned_lengths_match() {
        bail!("cannot merge alignments with inconsistent row lengths");
    }

    let left_profile = left.profile_tokens();
    let right_profile = right.profile_tokens();
    let trace = needleman_wunsch(&left_profile, &right_profile, scheme);

    let mut output_rows = vec![Vec::with_capacity(trace.len()); left.rows.len() + right.rows.len()];
    let mut left_column = 0usize;
    let mut right_column = 0usize;

    for op in trace {
        match op {
            TraceOp::Align => {
                for (row_index, row) in left.rows.iter().enumerate() {
                    output_rows[row_index].push(row[left_column]);
                }
                for (row_offset, row) in right.rows.iter().enumerate() {
                    output_rows[left.rows.len() + row_offset].push(row[right_column]);
                }
                left_column += 1;
                right_column += 1;
            }
            TraceOp::GapInLeft => {
                for row_index in 0..left.rows.len() {
                    output_rows[row_index].push(b'-');
                }
                for (row_offset, row) in right.rows.iter().enumerate() {
                    output_rows[left.rows.len() + row_offset].push(row[right_column]);
                }
                right_column += 1;
            }
            TraceOp::GapInRight => {
                for (row_index, row) in left.rows.iter().enumerate() {
                    output_rows[row_index].push(row[left_column]);
                }
                for row_index in left.rows.len()..output_rows.len() {
                    output_rows[row_index].push(b'-');
                }
                left_column += 1;
            }
        }
    }

    let mut ids = left.ids.clone();
    ids.extend(right.ids.clone());

    Ok(Alignment { ids, rows: output_rows })
}

fn needleman_wunsch(left: &[u8], right: &[u8], scheme: ScoreScheme) -> Vec<TraceOp> {
    let width = right.len() + 1;
    let mut scores = vec![0i32; (left.len() + 1) * width];
    let mut trace = vec![TraceOp::Align; (left.len() + 1) * width];

    for i in 1..=left.len() {
        scores[i * width] = scores[(i - 1) * width] + scheme.gap_penalty;
        trace[i * width] = TraceOp::GapInRight;
    }
    for j in 1..=right.len() {
        scores[j] = scores[j - 1] + scheme.gap_penalty;
        trace[j] = TraceOp::GapInLeft;
    }

    for i in 1..=left.len() {
        for j in 1..=right.len() {
            let diagonal = scores[(i - 1) * width + (j - 1)]
                + substitution_score(left[i - 1], right[j - 1], scheme);
            let gap_left = scores[i * width + (j - 1)] + scheme.gap_penalty;
            let gap_right = scores[(i - 1) * width + j] + scheme.gap_penalty;

            let (best_score, best_trace) = if diagonal >= gap_left && diagonal >= gap_right {
                (diagonal, TraceOp::Align)
            } else if gap_left >= gap_right {
                (gap_left, TraceOp::GapInLeft)
            } else {
                (gap_right, TraceOp::GapInRight)
            };

            scores[i * width + j] = best_score;
            trace[i * width + j] = best_trace;
        }
    }

    let mut ops = Vec::with_capacity(left.len() + right.len());
    let mut i = left.len();
    let mut j = right.len();

    while i > 0 || j > 0 {
        let op = trace[i * width + j];
        ops.push(op);
        match op {
            TraceOp::Align => {
                i -= 1;
                j -= 1;
            }
            TraceOp::GapInLeft => {
                j -= 1;
            }
            TraceOp::GapInRight => {
                i -= 1;
            }
        }
    }

    ops.reverse();
    ops
}

fn substitution_score(left: u8, right: u8, scheme: ScoreScheme) -> i32 {
    match (left, right) {
        (b'-', b'-') => 0,
        (b'-', _) | (_, b'-') => scheme.gap_penalty,
        _ if left == right => scheme.match_score,
        _ => scheme.mismatch_score,
    }
}

#[cfg(test)]
mod tests {
    use super::{polish_alignment, progressive_align, ScoreScheme};
    use crate::fasta::Sequence;
    use crate::guide::GuideTree;

    fn scheme() -> ScoreScheme {
        ScoreScheme {
            match_score: 3,
            mismatch_score: -1,
            gap_penalty: -2,
        }
    }

    #[test]
    fn progressive_alignment_keeps_rows_rectangular() {
        let sequences = vec![
            Sequence {
                id: "a".into(),
                residues: b"ACDEFG".to_vec(),
            },
            Sequence {
                id: "b".into(),
                residues: b"ACDFFG".to_vec(),
            },
            Sequence {
                id: "c".into(),
                residues: b"ACDFG".to_vec(),
            },
        ];
        let tree = GuideTree::Merge {
            left: Box::new(GuideTree::Leaf(0)),
            right: Box::new(GuideTree::Merge {
                left: Box::new(GuideTree::Leaf(1)),
                right: Box::new(GuideTree::Leaf(2)),
            }),
        };

        let alignment = progressive_align(&tree, &sequences, scheme()).unwrap();
        assert_eq!(alignment.rows.len(), 3);
        assert!(alignment.rows.iter().all(|row| row.len() == alignment.rows[0].len()));
    }

    #[test]
    fn polishing_preserves_all_sequences() {
        let sequences = vec![
            Sequence {
                id: "a".into(),
                residues: b"ACDEFG".to_vec(),
            },
            Sequence {
                id: "b".into(),
                residues: b"ACDEYG".to_vec(),
            },
            Sequence {
                id: "c".into(),
                residues: b"ACDFG".to_vec(),
            },
        ];
        let tree = GuideTree::Merge {
            left: Box::new(GuideTree::Leaf(0)),
            right: Box::new(GuideTree::Merge {
                left: Box::new(GuideTree::Leaf(1)),
                right: Box::new(GuideTree::Leaf(2)),
            }),
        };

        let alignment = progressive_align(&tree, &sequences, scheme()).unwrap();
        let polished = polish_alignment(alignment, &sequences, 1, scheme()).unwrap();

        assert_eq!(polished.ids.len(), 3);
        assert!(polished.ids.iter().any(|id| id == "a"));
        assert!(polished.ids.iter().any(|id| id == "b"));
        assert!(polished.ids.iter().any(|id| id == "c"));
    }
}
