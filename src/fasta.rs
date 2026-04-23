use std::collections::HashMap;
use std::fs;
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::aligner::Alignment;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Sequence {
    pub id: String,
    pub residues: Vec<u8>,
}

impl Sequence {
    pub fn len(&self) -> usize {
        self.residues.len()
    }
}

/// Dispatch to the correct reader based on the file extension.
/// Supported input formats:
///   FASTA   — .fa .fasta .faa .fna (default)
///   FASTQ   — .fastq .fq
///   Stockholm — .sto .stk .stockholm
pub fn read_sequences(path: &Path) -> Result<Vec<Sequence>> {
    match path.extension().and_then(|e| e.to_str()) {
        Some("fastq") | Some("fq") => read_fastq(path),
        Some("sto") | Some("stk") | Some("stockholm") => read_stockholm(path),
        _ => read_fasta(path),
    }
}

pub fn read_fasta(path: &Path) -> Result<Vec<Sequence>> {
    let raw = fs::read_to_string(path)
        .with_context(|| format!("failed to read FASTA input {}", path.display()))?;
    let mut sequences = Vec::new();
    let mut current_id: Option<String> = None;
    let mut current_residues = Vec::new();

    for line in raw.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix('>') {
            if let Some(id) = current_id.take() {
                if current_residues.is_empty() {
                    bail!("record {id} in {} had no residues", path.display());
                }
                sequences.push(Sequence {
                    id,
                    residues: std::mem::take(&mut current_residues),
                });
            }
            let id = rest.trim();
            if id.is_empty() {
                bail!("found a FASTA header without an identifier in {}", path.display());
            }
            current_id = Some(id.to_string());
            continue;
        }

        if current_id.is_none() {
            bail!(
                "encountered residues before the first FASTA header in {}",
                path.display()
            );
        }

        for byte in trimmed.bytes().filter(|b| !b.is_ascii_whitespace()) {
            if byte != b'-' {
                current_residues.push(byte.to_ascii_uppercase());
            }
        }
    }

    if let Some(id) = current_id {
        if current_residues.is_empty() {
            bail!("record {id} in {} had no residues", path.display());
        }
        sequences.push(Sequence {
            id,
            residues: current_residues,
        });
    }

    if sequences.is_empty() {
        bail!("{} did not contain any FASTA records", path.display());
    }

    Ok(sequences)
}

pub fn read_fastq(path: &Path) -> Result<Vec<Sequence>> {
    let raw = fs::read_to_string(path)
        .with_context(|| format!("failed to read FASTQ input {}", path.display()))?;

    let mut sequences = Vec::new();
    let mut lines = raw.lines().filter(|l| !l.trim().is_empty());

    loop {
        let header = match lines.next() {
            Some(h) => h.trim(),
            None => break,
        };

        let id = header
            .strip_prefix('@')
            .with_context(|| {
                format!(
                    "expected '@' header line in {}, got: {}",
                    path.display(),
                    header
                )
            })?
            .split_whitespace()
            .next()
            .unwrap_or("")
            .trim();

        if id.is_empty() {
            bail!("FASTQ record in {} has an empty identifier", path.display());
        }

        let seq_line = lines.next().with_context(|| {
            format!(
                "FASTQ record '{}' in {} is missing its sequence line",
                id,
                path.display()
            )
        })?;

        let residues: Vec<u8> = seq_line
            .trim()
            .bytes()
            .filter(|b| !b.is_ascii_whitespace() && *b != b'-')
            .map(|b| b.to_ascii_uppercase())
            .collect();

        if residues.is_empty() {
            bail!(
                "FASTQ record '{}' in {} has an empty sequence",
                id,
                path.display()
            );
        }

        let sep = lines.next().with_context(|| {
            format!(
                "FASTQ record '{}' in {} is missing the '+' separator line",
                id,
                path.display()
            )
        })?;
        if !sep.trim().starts_with('+') {
            bail!(
                "FASTQ record '{}' in {} expected '+' separator, got: {}",
                id,
                path.display(),
                sep.trim()
            );
        }

        // quality line — not used for alignment, just consumed
        lines.next().with_context(|| {
            format!(
                "FASTQ record '{}' in {} is missing its quality line",
                id,
                path.display()
            )
        })?;

        sequences.push(Sequence {
            id: id.to_string(),
            residues,
        });
    }

    if sequences.is_empty() {
        bail!("{} did not contain any FASTQ records", path.display());
    }

    Ok(sequences)
}

pub fn read_stockholm(path: &Path) -> Result<Vec<Sequence>> {
    let raw = fs::read_to_string(path)
        .with_context(|| format!("failed to read Stockholm input {}", path.display()))?;

    let mut order: Vec<String> = Vec::new();
    let mut map: HashMap<String, Vec<u8>> = HashMap::new();

    for line in raw.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') || trimmed == "//" {
            continue;
        }

        let mut parts = trimmed.splitn(2, char::is_whitespace);
        let id = parts.next().unwrap_or("").trim();
        let seq = parts.next().unwrap_or("").trim();

        if id.is_empty() || seq.is_empty() {
            continue;
        }

        let residues: Vec<u8> = seq
            .bytes()
            .filter(|b| !b.is_ascii_whitespace() && *b != b'.' && *b != b'-')
            .map(|b| b.to_ascii_uppercase())
            .collect();

        if !map.contains_key(id) {
            order.push(id.to_string());
            map.insert(id.to_string(), Vec::new());
        }
        map.get_mut(id).unwrap().extend_from_slice(&residues);
    }

    let sequences: Vec<Sequence> = order
        .into_iter()
        .filter_map(|id| {
            let residues = map.remove(&id).unwrap_or_default();
            if residues.is_empty() {
                None
            } else {
                Some(Sequence { id, residues })
            }
        })
        .collect();

    if sequences.is_empty() {
        bail!("{} did not contain any Stockholm sequence records", path.display());
    }

    Ok(sequences)
}

/// Write an alignment. Output format is chosen by the file extension:
///   Clustal — .aln .clustal
///   FASTA   — everything else (default)
pub fn write_alignment(path: &Path, alignment: &Alignment) -> Result<()> {
    if alignment.rows.is_empty() {
        bail!("refusing to write an empty alignment to {}", path.display());
    }

    match path.extension().and_then(|e| e.to_str()) {
        Some("aln") | Some("clustal") => write_clustal(path, alignment),
        _ => write_fasta(path, alignment),
    }
}

fn write_fasta(path: &Path, alignment: &Alignment) -> Result<()> {
    let mut output = String::new();
    for (id, row) in alignment.ids.iter().zip(&alignment.rows) {
        output.push('>');
        output.push_str(id);
        output.push('\n');
        for chunk in row.chunks(80) {
            output.push_str(
                std::str::from_utf8(chunk).context("alignment contained non-UTF8 bytes")?,
            );
            output.push('\n');
        }
    }

    fs::write(path, output)
        .with_context(|| format!("failed to write aligned FASTA to {}", path.display()))?;
    Ok(())
}

fn write_clustal(path: &Path, alignment: &Alignment) -> Result<()> {
    const BLOCK: usize = 60;

    let id_width = alignment
        .ids
        .iter()
        .map(|id| id.len())
        .max()
        .unwrap_or(0)
        + 4;

    let columns = alignment.column_count();
    let mut output = String::from("CLUSTAL W (gpumsa) multiple sequence alignment\n\n");

    let mut col = 0;
    while col < columns {
        let end = (col + BLOCK).min(columns);
        for (id, row) in alignment.ids.iter().zip(&alignment.rows) {
            let slice = std::str::from_utf8(&row[col..end])
                .context("alignment contained non-UTF8 bytes")?;
            output.push_str(&format!("{:<width$}{}\n", id, slice, width = id_width));
        }
        output.push('\n');
        col = end;
    }

    fs::write(path, output)
        .with_context(|| format!("failed to write Clustal alignment to {}", path.display()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{read_fasta, read_fastq, read_stockholm, write_alignment};
    use crate::aligner::Alignment;
    use std::fs;

    fn tmp(tag: &str, ext: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("gpumsa-{}-{}.{}", tag, std::process::id(), ext))
    }

    // ── FASTA happy-path ─────────────────────────────────────────────────────

    #[test]
    fn parses_multiline_fasta() {
        let path = tmp("fasta-parse", "fa");
        fs::write(&path, ">seq1\nacde\nfg\n>seq2\nMK-T\n").unwrap();

        let records = read_fasta(&path).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].residues, b"ACDEFG");
        assert_eq!(records[1].residues, b"MKT");

        let _ = fs::remove_file(path);
    }

    #[test]
    fn fasta_strips_gaps_from_prealigned_input() {
        let path = tmp("fasta-gaps", "fa");
        fs::write(&path, ">s1\nMK--TA\n>s2\nMK--TA\n").unwrap();
        let records = read_fasta(&path).unwrap();
        assert_eq!(records[0].residues, b"MKTA");
        assert_eq!(records[1].residues, b"MKTA");
        let _ = fs::remove_file(path);
    }

    #[test]
    fn fasta_preserves_ambiguous_amino_acids() {
        let path = tmp("fasta-ambig", "fa");
        fs::write(&path, ">s1\nxbzuo\n").unwrap();
        let records = read_fasta(&path).unwrap();
        assert_eq!(records[0].residues, b"XBZUO");
        let _ = fs::remove_file(path);
    }

    #[test]
    fn fasta_ignores_blank_lines_between_records() {
        let path = tmp("fasta-blanks", "fa");
        fs::write(&path, "\n>s1\nMKT\n\n>s2\nACG\n\n").unwrap();
        let records = read_fasta(&path).unwrap();
        assert_eq!(records.len(), 2);
        let _ = fs::remove_file(path);
    }

    // ── FASTA error cases ────────────────────────────────────────────────────

    #[test]
    fn fasta_error_on_empty_file() {
        let path = tmp("fasta-empty", "fa");
        fs::write(&path, "").unwrap();
        assert!(read_fasta(&path).is_err());
        let _ = fs::remove_file(path);
    }

    #[test]
    fn fasta_error_on_whitespace_only_file() {
        let path = tmp("fasta-ws", "fa");
        fs::write(&path, "\n   \n\t\n").unwrap();
        assert!(read_fasta(&path).is_err());
        let _ = fs::remove_file(path);
    }

    #[test]
    fn fasta_error_residues_before_first_header() {
        let path = tmp("fasta-nohdr", "fa");
        fs::write(&path, "ACGT\n>seq1\nACGT\n").unwrap();
        assert!(read_fasta(&path).is_err());
        let _ = fs::remove_file(path);
    }

    #[test]
    fn fasta_error_empty_header_id() {
        let path = tmp("fasta-emptyid", "fa");
        fs::write(&path, ">\nACGT\n").unwrap();
        assert!(read_fasta(&path).is_err());
        let _ = fs::remove_file(path);
    }

    #[test]
    fn fasta_error_record_with_no_residues() {
        let path = tmp("fasta-nores", "fa");
        fs::write(&path, ">s1\n>s2\nACGT\n").unwrap();
        assert!(read_fasta(&path).is_err());
        let _ = fs::remove_file(path);
    }

    #[test]
    fn fasta_error_last_record_with_no_residues() {
        let path = tmp("fasta-lastempty", "fa");
        fs::write(&path, ">s1\nACGT\n>s2\n").unwrap();
        assert!(read_fasta(&path).is_err());
        let _ = fs::remove_file(path);
    }

    // ── FASTQ happy-path ─────────────────────────────────────────────────────

    #[test]
    fn parses_fastq() {
        let path = tmp("fastq-parse", "fastq");
        fs::write(
            &path,
            "@read1 some description\nMKTAYIAK\n+\nIIIIIIII\n@read2\nACDEFG\n+read2\n!!!!!!\n",
        )
        .unwrap();

        let records = read_fastq(&path).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "read1");
        assert_eq!(records[0].residues, b"MKTAYIAK");
        assert_eq!(records[1].id, "read2");
        assert_eq!(records[1].residues, b"ACDEFG");

        let _ = fs::remove_file(path);
    }

    #[test]
    fn fastq_strips_gaps_and_uppercases() {
        let path = tmp("fastq-upper", "fq");
        fs::write(&path, "@s1\nakt-gcy\n+\n!!!!!!!\n").unwrap();
        let records = read_fastq(&path).unwrap();
        assert_eq!(records[0].residues, b"AKTGCY");
        let _ = fs::remove_file(path);
    }

    #[test]
    fn fastq_multiword_header_uses_first_word_as_id() {
        let path = tmp("fastq-desc", "fq");
        fs::write(&path, "@read1 extra description info\nMKT\n+\n!!!\n").unwrap();
        let records = read_fastq(&path).unwrap();
        assert_eq!(records[0].id, "read1");
        let _ = fs::remove_file(path);
    }

    // ── FASTQ error cases ────────────────────────────────────────────────────

    #[test]
    fn fastq_error_on_empty_file() {
        let path = tmp("fastq-empty", "fastq");
        fs::write(&path, "").unwrap();
        assert!(read_fastq(&path).is_err());
        let _ = fs::remove_file(path);
    }

    #[test]
    fn fastq_error_missing_at_prefix() {
        let path = tmp("fastq-noat", "fastq");
        fs::write(&path, "read1\nMKT\n+\n!!!\n").unwrap();
        assert!(read_fastq(&path).is_err());
        let _ = fs::remove_file(path);
    }

    #[test]
    fn fastq_error_missing_plus_separator() {
        let path = tmp("fastq-noplus", "fastq");
        fs::write(&path, "@read1\nMKT\nMKT\n!!!\n").unwrap();
        assert!(read_fastq(&path).is_err());
        let _ = fs::remove_file(path);
    }

    #[test]
    fn fastq_error_truncated_after_sequence() {
        let path = tmp("fastq-trunc", "fastq");
        fs::write(&path, "@read1\nMKT\n").unwrap();
        assert!(read_fastq(&path).is_err());
        let _ = fs::remove_file(path);
    }

    #[test]
    fn fastq_error_truncated_after_plus() {
        let path = tmp("fastq-trunc2", "fastq");
        fs::write(&path, "@read1\nMKT\n+\n").unwrap();
        assert!(read_fastq(&path).is_err());
        let _ = fs::remove_file(path);
    }

    // ── Stockholm happy-path ─────────────────────────────────────────────────

    #[test]
    fn parses_stockholm() {
        let path = tmp("sto-parse", "sto");
        fs::write(
            &path,
            "# STOCKHOLM 1.0\n#=GF ID  test\nseq1  MKTA-YIAK\nseq2  MKTG.YIAK\n//\n",
        )
        .unwrap();

        let records = read_stockholm(&path).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].residues, b"MKTAYIAK");
        assert_eq!(records[1].id, "seq2");
        assert_eq!(records[1].residues, b"MKTGYIAK");

        let _ = fs::remove_file(path);
    }

    #[test]
    fn stockholm_interleaved_multi_block() {
        let path = tmp("sto-interleaved", "sto");
        // Two alignment blocks — residues should be concatenated per sequence
        fs::write(
            &path,
            "# STOCKHOLM 1.0\nseq1  MKTA\nseq2  MKTG\n\nseq1  YIAK\nseq2  YIAK\n//\n",
        )
        .unwrap();
        let records = read_stockholm(&path).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].residues, b"MKTAYIAK");
        assert_eq!(records[1].residues, b"MKTGYIAK");
        let _ = fs::remove_file(path);
    }

    #[test]
    fn stockholm_skips_markup_and_comment_lines() {
        let path = tmp("sto-markup", "sto");
        fs::write(
            &path,
            "# STOCKHOLM 1.0\n#=GF AC PF00001\n#=GC SS_cons ....\nseq1  MKTA\n//\n",
        )
        .unwrap();
        let records = read_stockholm(&path).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id, "seq1");
        let _ = fs::remove_file(path);
    }

    // ── Stockholm error cases ────────────────────────────────────────────────

    #[test]
    fn stockholm_error_on_empty_file() {
        let path = tmp("sto-empty", "sto");
        fs::write(&path, "").unwrap();
        assert!(read_stockholm(&path).is_err());
        let _ = fs::remove_file(path);
    }

    #[test]
    fn stockholm_error_only_markup_no_sequences() {
        let path = tmp("sto-nodata", "sto");
        fs::write(&path, "# STOCKHOLM 1.0\n#=GF ID  empty\n//\n").unwrap();
        assert!(read_stockholm(&path).is_err());
        let _ = fs::remove_file(path);
    }

    #[test]
    fn stockholm_error_sequences_that_are_only_gaps() {
        let path = tmp("sto-allgap", "sto");
        fs::write(&path, "# STOCKHOLM 1.0\nseq1  ----\nseq2  ....\n//\n").unwrap();
        assert!(read_stockholm(&path).is_err());
        let _ = fs::remove_file(path);
    }

    // ── Output format tests ──────────────────────────────────────────────────

    #[test]
    fn fasta_output_wraps_at_80_chars() {
        let path = tmp("fasta-wrap", "fa");
        let long_row: Vec<u8> = b"A".repeat(200).to_vec();
        let alignment = Alignment {
            ids: vec!["s1".into()],
            rows: vec![long_row],
        };
        write_alignment(&path, &alignment).unwrap();
        let content = fs::read_to_string(&path).unwrap();
        for line in content.lines() {
            if !line.starts_with('>') {
                assert!(line.len() <= 80, "line too long: {} chars", line.len());
            }
        }
        let _ = fs::remove_file(path);
    }

    #[test]
    fn clustal_output_splits_into_60_char_blocks() {
        let path = tmp("clustal-blocks", "aln");
        let row: Vec<u8> = b"M".repeat(130).to_vec();
        let alignment = Alignment {
            ids: vec!["seq1".into(), "seq2".into()],
            rows: vec![row.clone(), row],
        };
        write_alignment(&path, &alignment).unwrap();
        let content = fs::read_to_string(&path).unwrap();
        // Every non-header, non-blank data line should have at most 60 residue chars
        // (plus the padded ID prefix)
        for line in content.lines() {
            if line.is_empty() || line.starts_with("CLUSTAL") {
                continue;
            }
            let residues = line.trim_start_matches(|c: char| !c.is_ascii_uppercase() && c != '-');
            assert!(residues.len() <= 60, "block line too wide: {}", residues.len());
        }
        let _ = fs::remove_file(path);
    }

    #[test]
    fn write_alignment_errors_on_empty_alignment() {
        let path = tmp("empty-aln", "fa");
        let alignment = Alignment { ids: vec![], rows: vec![] };
        assert!(write_alignment(&path, &alignment).is_err());
        let _ = fs::remove_file(path);
    }
}
