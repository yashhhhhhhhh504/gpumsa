use std::path::PathBuf;

use anyhow::Result;
use clap::{Parser, ValueEnum};
use gpumsa::{detect_gpu, run, BackendKind, Config, ScoreScheme};

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
enum BackendChoice {
    Auto,
    Cpu,
    Gpu,
}

impl From<BackendChoice> for BackendKind {
    fn from(value: BackendChoice) -> Self {
        match value {
            BackendChoice::Auto => BackendKind::Auto,
            BackendChoice::Cpu => BackendKind::Cpu,
            BackendChoice::Gpu => BackendKind::Gpu,
        }
    }
}

#[derive(Debug, Parser)]
#[command(
    name = "gpumsa",
    about = "Prototype GPU-accelerated multiple sequence aligner for protein families"
)]
struct Cli {
    /// Probe the system for a compatible GPU and print adapter info, then exit
    #[arg(long)]
    check_gpu: bool,

    /// Input sequence file (.fa/.fasta/.faa/.fna, .fastq/.fq, .sto/.stk)
    #[arg(short, long, value_name = "FILE")]
    input: Option<PathBuf>,
    /// Output alignment file (.fa/.fasta for FASTA, .aln/.clustal for Clustal)
    #[arg(short, long, value_name = "FILE")]
    output: Option<PathBuf>,
    #[arg(long, value_enum, default_value = "auto")]
    backend: BackendChoice,
    #[arg(long, default_value_t = 1)]
    refine_passes: usize,
    #[arg(long, default_value_t = 4)]
    match_score: i32,
    #[arg(long, default_value_t = -2)]
    mismatch_score: i32,
    #[arg(long, default_value_t = -6)]
    gap_penalty: i32,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.check_gpu {
        println!("Probing system for GPU adapters …\n");
        match detect_gpu() {
            Some(info) => {
                println!("GPU detected:");
                println!("  Name        : {}", info.name);
                println!("  Backend     : {}", info.backend);
                println!("  Device type : {}", info.device_type);
                println!("  Driver      : {}", info.driver);
                println!("\nGPU is available. Use --backend gpu or --backend auto to accelerate alignments.");
            }
            None => {
                println!("No compatible GPU adapter found.");
                println!("The tool will use CPU-only mode (--backend cpu).");
            }
        }
        return Ok(());
    }

    let input = cli
        .input
        .ok_or_else(|| anyhow::anyhow!("--input is required for alignment (use --check-gpu to probe GPU)"))?;
    let output = cli
        .output
        .ok_or_else(|| anyhow::anyhow!("--output is required for alignment"))?;

    let summary = run(Config {
        input,
        output,
        backend: cli.backend.into(),
        refine_passes: cli.refine_passes,
        score_scheme: ScoreScheme {
            match_score: cli.match_score,
            mismatch_score: cli.mismatch_score,
            gap_penalty: cli.gap_penalty,
        },
    })?;

    println!(
        "Aligned {} sequences into {} columns using the {} similarity backend. Wrote {}.",
        summary.sequences,
        summary.columns,
        summary.backend,
        summary.output.display()
    );

    Ok(())
}
