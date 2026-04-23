pub mod aligner;
pub mod fasta;
pub mod gpu;
pub mod guide;
pub mod pipeline;
pub mod similarity;

// Core pipeline
pub use pipeline::{run, Config, RunSummary};

// Configuration types
pub use aligner::ScoreScheme;
pub use similarity::BackendKind;

// GPU probe
pub use gpu::{detect_gpu, GpuInfo};

// Format I/O — exposed so library consumers can load sequences or write
// alignments without going through `run()`.
pub use fasta::{read_fasta, read_fastq, read_sequences, read_stockholm, write_alignment, Sequence};
