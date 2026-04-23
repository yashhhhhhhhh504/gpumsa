/// End-to-end library usage demo.
///
/// Run with:
///   cargo run --example align_demo
///
/// This example exercises the full public API:
///   - `detect_gpu` / `GpuInfo`
///   - `read_sequences` (auto-dispatch), `read_fasta`, `read_fastq`, `read_stockholm`
///   - `run` + `Config` + `ScoreScheme` + `BackendKind`
///   - `write_alignment` (FASTA and Clustal W output)
///   - `Sequence` (direct construction)
use std::path::Path;

use anyhow::Result;
use gpumsa::{
    detect_gpu, read_fasta, read_fastq, read_sequences, read_stockholm, run, write_alignment,
    BackendKind, Config, Sequence, ScoreScheme,
};

fn main() -> Result<()> {
    // ── 1. GPU probe ─────────────────────────────────────────────────────────
    println!("=== GPU probe ===");
    match detect_gpu() {
        Some(info) => println!("  GPU: {} ({}, {})", info.name, info.backend, info.device_type),
        None => println!("  No GPU found — CPU fallback will be used"),
    }

    // ── 2. Read sequences via the auto-dispatch function ──────────────────────
    println!("\n=== Auto-dispatch reader ===");
    let fasta_seqs = read_sequences(Path::new("examples/toy_proteins.fasta"))?;
    println!("  FASTA:     {} sequences", fasta_seqs.len());

    let fastq_seqs = read_sequences(Path::new("examples/toy_proteins.fastq"))?;
    println!("  FASTQ:     {} sequences", fastq_seqs.len());

    let sto_seqs = read_sequences(Path::new("examples/toy_proteins.sto"))?;
    println!("  Stockholm: {} sequences", sto_seqs.len());

    // All three files contain the same sequences — verify residues match
    assert_eq!(fasta_seqs.len(), fastq_seqs.len());
    assert_eq!(fasta_seqs.len(), sto_seqs.len());
    for (fa, fq) in fasta_seqs.iter().zip(&fastq_seqs) {
        assert_eq!(fa.residues, fq.residues, "FASTA/FASTQ mismatch for {}", fa.id);
    }
    for (fa, sto) in fasta_seqs.iter().zip(&sto_seqs) {
        assert_eq!(fa.residues, sto.residues, "FASTA/Stockholm mismatch for {}", fa.id);
    }
    println!("  All three files produce identical residues ✓");

    // ── 3. Use specific readers directly ──────────────────────────────────────
    println!("\n=== Format-specific readers ===");
    let via_fasta = read_fasta(Path::new("examples/toy_proteins.fasta"))?;
    let via_fastq = read_fastq(Path::new("examples/toy_proteins.fastq"))?;
    let via_sto = read_stockholm(Path::new("examples/toy_proteins.sto"))?;
    println!("  read_fasta:     {} seqs, first id = {}", via_fasta.len(), via_fasta[0].id);
    println!("  read_fastq:     {} seqs, first id = {}", via_fastq.len(), via_fastq[0].id);
    println!("  read_stockholm: {} seqs, first id = {}", via_sto.len(), via_sto[0].id);

    // ── 4. Construct Sequence directly ────────────────────────────────────────
    println!("\n=== Direct Sequence construction ===");
    let custom = Sequence {
        id: "custom_kinase".into(),
        residues: b"MKTAYIAKQRQISFVK".to_vec(),
    };
    println!("  Sequence {{ id: {}, len: {} }}", custom.id, custom.len());

    // ── 5. Run alignment (FASTA output) ───────────────────────────────────────
    println!("\n=== Alignment: FASTA → FASTA ===");
    let fasta_out = std::env::temp_dir().join("demo_aligned.fasta");
    let scheme = ScoreScheme { match_score: 4, mismatch_score: -2, gap_penalty: -6 };

    let summary = run(Config {
        input: "examples/toy_proteins.fasta".into(),
        output: fasta_out.clone(),
        backend: BackendKind::Auto,
        refine_passes: 1,
        score_scheme: scheme,
    })?;

    println!(
        "  Aligned {} sequences into {} columns via {} backend",
        summary.sequences, summary.columns, summary.backend
    );
    println!("  Output: {}", summary.output.display());

    // ── 6. Run alignment (Clustal output) ────────────────────────────────────
    println!("\n=== Alignment: FASTQ → Clustal W ===");
    let clustal_out = std::env::temp_dir().join("demo_aligned.aln");

    let summary2 = run(Config {
        input: "examples/toy_proteins.fastq".into(),
        output: clustal_out.clone(),
        backend: BackendKind::Cpu,
        refine_passes: 1,
        score_scheme: scheme,
    })?;

    println!(
        "  Aligned {} sequences into {} columns via {} backend",
        summary2.sequences, summary2.columns, summary2.backend
    );
    let clustal_content = std::fs::read_to_string(&clustal_out)?;
    assert!(clustal_content.starts_with("CLUSTAL W"), "Clustal header missing");
    println!("  Clustal header present ✓");

    // ── 7. write_alignment directly ──────────────────────────────────────────
    println!("\n=== Direct write_alignment ===");
    // Re-read the FASTA output we just wrote and write it as Clustal W
    let realigned = read_fasta(&fasta_out)?;
    // Build an Alignment from the re-read sequences (no re-alignment, just reformat)
    use gpumsa::aligner::Alignment;
    let manual_alignment = Alignment {
        ids: realigned.iter().map(|s| s.id.clone()).collect(),
        rows: realigned.iter().map(|s| s.residues.clone()).collect(),
    };
    let reformat_out = std::env::temp_dir().join("demo_reformatted.aln");
    write_alignment(&reformat_out, &manual_alignment)?;
    let reformat_content = std::fs::read_to_string(&reformat_out)?;
    assert!(reformat_content.starts_with("CLUSTAL W"));
    println!("  Reformatted FASTA → Clustal W via write_alignment ✓");

    // ── 8. Stockholm alignment ────────────────────────────────────────────────
    println!("\n=== Alignment: Stockholm → FASTA ===");
    let sto_out = std::env::temp_dir().join("demo_from_sto.fasta");
    let summary3 = run(Config {
        input: "examples/toy_proteins.sto".into(),
        output: sto_out.clone(),
        backend: BackendKind::Cpu,
        refine_passes: 1,
        score_scheme: scheme,
    })?;
    println!(
        "  Aligned {} sequences into {} columns",
        summary3.sequences, summary3.columns
    );
    assert_eq!(summary3.sequences, 4);

    println!("\nAll checks passed.");
    Ok(())
}
