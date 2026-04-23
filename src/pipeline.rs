use std::path::PathBuf;

use anyhow::{bail, Result};

use crate::aligner::{polish_alignment, progressive_align, ScoreScheme};
use crate::fasta::{read_sequences, write_alignment};
use crate::guide::build_guide_tree;
use crate::similarity::{compute_similarity_matrix, BackendKind};

#[derive(Clone, Debug)]
pub struct Config {
    pub input: PathBuf,
    pub output: PathBuf,
    pub backend: BackendKind,
    pub refine_passes: usize,
    pub score_scheme: ScoreScheme,
}

#[derive(Clone, Debug)]
pub struct RunSummary {
    pub sequences: usize,
    pub columns: usize,
    pub backend: &'static str,
    pub output: PathBuf,
}

pub fn run(config: Config) -> Result<RunSummary> {
    let sequences = read_sequences(&config.input)?;
    if sequences.is_empty() {
        bail!("input FASTA contained no sequences");
    }

    let scored = compute_similarity_matrix(&sequences, config.backend)?;
    let tree = build_guide_tree(&scored.matrix, sequences.len());
    let aligned = progressive_align(&tree, &sequences, config.score_scheme)?;
    let polished = polish_alignment(aligned, &sequences, config.refine_passes, config.score_scheme)?;
    write_alignment(&config.output, &polished)?;

    Ok(RunSummary {
        sequences: sequences.len(),
        columns: polished.column_count(),
        backend: scored.backend,
        output: config.output,
    })
}

#[cfg(test)]
mod tests {
    use super::{run, Config};
    use crate::aligner::ScoreScheme;
    use crate::similarity::BackendKind;
    use std::fs;

    fn default_scheme() -> ScoreScheme {
        ScoreScheme {
            match_score: 4,
            mismatch_score: -2,
            gap_penalty: -6,
        }
    }

    fn base_path(tag: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("gpumsa-{}-{}", tag, std::process::id()))
    }

    #[test]
    fn fasta_input_fasta_output() {
        let input = base_path("fasta-in").with_extension("fa");
        let output = base_path("fasta-out").with_extension("fa");
        fs::write(
            &input,
            ">seq1\nMKTAYIAKQRQISFVKSHFSRQDILDLWQ\n>seq2\nMKTAYIAKQRQISFVKSHFSRNDILDLWQ\n>seq3\nMKTTYIAKQRQISFVKAHFSRQDILDLWQ\n",
        )
        .unwrap();

        let summary = run(Config {
            input: input.clone(),
            output: output.clone(),
            backend: BackendKind::Cpu,
            refine_passes: 1,
            score_scheme: default_scheme(),
        })
        .unwrap();

        let aligned = fs::read_to_string(&output).unwrap();
        assert_eq!(summary.sequences, 3);
        assert!(summary.columns >= 29);
        assert!(aligned.contains(">seq1"));
        assert!(aligned.contains(">seq2"));
        assert!(aligned.contains(">seq3"));

        let _ = fs::remove_file(input);
        let _ = fs::remove_file(output);
    }

    #[test]
    fn fastq_input_fasta_output() {
        let input = base_path("fastq-in").with_extension("fastq");
        let output = base_path("fastq-out").with_extension("fa");
        fs::write(
            &input,
            "@seq1\nMKTAYIAKQRQISFVKSHFSRQDILDLWQ\n+\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\
             @seq2\nMKTAYIAKQRQISFVKSHFSRNDILDLWQ\n+\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\
             @seq3\nMKTTYIAKQRQISFVKAHFSRQDILDLWQ\n+\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
        )
        .unwrap();

        let summary = run(Config {
            input: input.clone(),
            output: output.clone(),
            backend: BackendKind::Cpu,
            refine_passes: 1,
            score_scheme: default_scheme(),
        })
        .unwrap();

        let aligned = fs::read_to_string(&output).unwrap();
        assert_eq!(summary.sequences, 3);
        assert!(summary.columns >= 29);
        assert!(aligned.contains(">seq1"));
        assert!(aligned.contains(">seq2"));
        assert!(aligned.contains(">seq3"));

        let _ = fs::remove_file(input);
        let _ = fs::remove_file(output);
    }

    #[test]
    fn stockholm_input_fasta_output() {
        let input = base_path("sto-in").with_extension("sto");
        let output = base_path("sto-out").with_extension("fa");
        fs::write(
            &input,
            "# STOCKHOLM 1.0\nseq1  MKTAYIAKQRQISFVKSHFSRQDILDLWQ\nseq2  MKTAYIAKQRQISFVKSHFSRNDILDLWQ\nseq3  MKTTYIAKQRQISFVKAHFSRQDILDLWQ\n//\n",
        )
        .unwrap();

        let summary = run(Config {
            input: input.clone(),
            output: output.clone(),
            backend: BackendKind::Cpu,
            refine_passes: 1,
            score_scheme: default_scheme(),
        })
        .unwrap();

        let aligned = fs::read_to_string(&output).unwrap();
        assert_eq!(summary.sequences, 3);
        assert!(summary.columns >= 29);
        assert!(aligned.contains(">seq1"));
        assert!(aligned.contains(">seq2"));
        assert!(aligned.contains(">seq3"));

        let _ = fs::remove_file(input);
        let _ = fs::remove_file(output);
    }

    #[test]
    fn fasta_input_clustal_output() {
        let input = base_path("clustal-in").with_extension("fa");
        let output = base_path("clustal-out").with_extension("aln");
        fs::write(
            &input,
            ">seq1\nMKTAYIAKQRQISFVKSHFSRQDILDLWQ\n>seq2\nMKTAYIAKQRQISFVKSHFSRNDILDLWQ\n>seq3\nMKTTYIAKQRQISFVKAHFSRQDILDLWQ\n",
        )
        .unwrap();

        run(Config {
            input: input.clone(),
            output: output.clone(),
            backend: BackendKind::Cpu,
            refine_passes: 1,
            score_scheme: default_scheme(),
        })
        .unwrap();

        let aligned = fs::read_to_string(&output).unwrap();
        assert!(aligned.starts_with("CLUSTAL W"));
        assert!(aligned.contains("seq1"));
        assert!(aligned.contains("seq2"));
        assert!(aligned.contains("seq3"));

        let _ = fs::remove_file(input);
        let _ = fs::remove_file(output);
    }

    #[test]
    fn single_sequence_passes_through_unchanged() {
        let input = base_path("single-in").with_extension("fa");
        let output = base_path("single-out").with_extension("fa");
        fs::write(&input, ">only\nMKTAYIAKQRQISFVK\n").unwrap();

        let summary = run(Config {
            input: input.clone(),
            output: output.clone(),
            backend: BackendKind::Cpu,
            refine_passes: 1,
            score_scheme: default_scheme(),
        })
        .unwrap();

        let aligned = fs::read_to_string(&output).unwrap();
        assert_eq!(summary.sequences, 1);
        assert_eq!(summary.columns, 16);
        assert!(aligned.contains("MKTAYIAKQRQISFVK"));

        let _ = fs::remove_file(input);
        let _ = fs::remove_file(output);
    }

    #[test]
    fn identical_sequences_align_without_gaps() {
        let seq = "MKTAYIAKQRQISFVK";
        let input = base_path("identical-in").with_extension("fa");
        let output = base_path("identical-out").with_extension("fa");
        fs::write(&input, format!(">s1\n{seq}\n>s2\n{seq}\n")).unwrap();

        let summary = run(Config {
            input: input.clone(),
            output: output.clone(),
            backend: BackendKind::Cpu,
            refine_passes: 1,
            score_scheme: default_scheme(),
        })
        .unwrap();

        assert_eq!(summary.columns, seq.len());
        let aligned = fs::read_to_string(&output).unwrap();
        assert!(!aligned.contains('-'));

        let _ = fs::remove_file(input);
        let _ = fs::remove_file(output);
    }

    #[test]
    fn unequal_length_sequences_shorter_gets_gaps() {
        let input = base_path("unequal-in").with_extension("fa");
        let output = base_path("unequal-out").with_extension("fa");
        fs::write(&input, ">short\nMK\n>long\nMKTAYIAKQRQISFVK\n").unwrap();

        let summary = run(Config {
            input: input.clone(),
            output: output.clone(),
            backend: BackendKind::Cpu,
            refine_passes: 1,
            score_scheme: default_scheme(),
        })
        .unwrap();

        assert!(summary.columns >= 16);
        let aligned = fs::read_to_string(&output).unwrap();
        assert!(aligned.contains('-'));

        let _ = fs::remove_file(input);
        let _ = fs::remove_file(output);
    }

    #[test]
    fn long_sequences_wrap_at_80_in_fasta_output() {
        let seq = "MKTAYIAKQRQISFVK".repeat(6); // 96 residues
        let input = base_path("long-in").with_extension("fa");
        let output = base_path("long-out").with_extension("fa");
        fs::write(&input, format!(">s1\n{seq}\n>s2\n{seq}\n")).unwrap();

        run(Config {
            input: input.clone(),
            output: output.clone(),
            backend: BackendKind::Cpu,
            refine_passes: 0,
            score_scheme: default_scheme(),
        })
        .unwrap();

        let raw = fs::read_to_string(&output).unwrap();
        for line in raw.lines() {
            if !line.starts_with('>') {
                assert!(line.len() <= 80, "residue line exceeds 80 chars: {}", line.len());
            }
        }

        let _ = fs::remove_file(input);
        let _ = fs::remove_file(output);
    }

    #[test]
    fn long_sequences_split_into_blocks_in_clustal_output() {
        let seq = "MKTAYIAKQRQISFVK".repeat(5); // 80 residues
        let input = base_path("long-clustal-in").with_extension("fa");
        let output = base_path("long-clustal-out").with_extension("aln");
        fs::write(&input, format!(">s1\n{seq}\n>s2\n{seq}\n")).unwrap();

        run(Config {
            input: input.clone(),
            output: output.clone(),
            backend: BackendKind::Cpu,
            refine_passes: 0,
            score_scheme: default_scheme(),
        })
        .unwrap();

        let raw = fs::read_to_string(&output).unwrap();
        assert!(raw.starts_with("CLUSTAL W"));
        // 80 residue columns at 60 per block → at least 2 data blocks
        let data_blocks = raw.split("\n\n").filter(|b| b.contains("s1")).count();
        assert!(data_blocks >= 2, "expected ≥2 Clustal blocks, got {}", data_blocks);

        let _ = fs::remove_file(input);
        let _ = fs::remove_file(output);
    }

    #[test]
    fn ambiguous_amino_acids_pass_through_pipeline() {
        let input = base_path("ambig-in").with_extension("fa");
        let output = base_path("ambig-out").with_extension("fa");
        // X B Z U are valid ambiguous amino acid codes
        fs::write(&input, ">s1\nMKTXBZU\n>s2\nMKTXBZU\n").unwrap();

        let summary = run(Config {
            input: input.clone(),
            output: output.clone(),
            backend: BackendKind::Cpu,
            refine_passes: 1,
            score_scheme: default_scheme(),
        })
        .unwrap();

        assert_eq!(summary.sequences, 2);
        let aligned = fs::read_to_string(&output).unwrap();
        assert!(aligned.contains("MKTXBZU"));

        let _ = fs::remove_file(input);
        let _ = fs::remove_file(output);
    }

    #[test]
    fn fastq_and_fasta_inputs_produce_same_alignment() {
        let fa_input = base_path("same-fa").with_extension("fa");
        let fq_input = base_path("same-fq").with_extension("fq");
        let fa_output = base_path("same-fa-out").with_extension("fa");
        let fq_output = base_path("same-fq-out").with_extension("fa");

        let seqs = ">seq1\nMKTAYIAKQRQISFVK\n>seq2\nMKTAYIAKQRQISFIK\n>seq3\nMKTTYIAKQRQISFVK\n";
        fs::write(&fa_input, seqs).unwrap();
        fs::write(
            &fq_input,
            "@seq1\nMKTAYIAKQRQISFVK\n+\n~~~~~~~~~~~~~~~~\n\
             @seq2\nMKTAYIAKQRQISFIK\n+\n~~~~~~~~~~~~~~~~\n\
             @seq3\nMKTTYIAKQRQISFVK\n+\n~~~~~~~~~~~~~~~~\n",
        )
        .unwrap();

        let cfg = |input: std::path::PathBuf, output: std::path::PathBuf| Config {
            input,
            output,
            backend: BackendKind::Cpu,
            refine_passes: 1,
            score_scheme: default_scheme(),
        };

        run(cfg(fa_input.clone(), fa_output.clone())).unwrap();
        run(cfg(fq_input.clone(), fq_output.clone())).unwrap();

        assert_eq!(
            fs::read_to_string(&fa_output).unwrap(),
            fs::read_to_string(&fq_output).unwrap()
        );

        let _ = fs::remove_file(fa_input);
        let _ = fs::remove_file(fq_input);
        let _ = fs::remove_file(fa_output);
        let _ = fs::remove_file(fq_output);
    }
}
