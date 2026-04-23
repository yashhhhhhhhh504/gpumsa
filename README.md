# gpumsa

> **GPU-accelerated multiple sequence aligner for large protein families**

`gpumsa` is a Rust library and CLI tool that uses [`wgpu`](https://wgpu.rs) compute shaders to accelerate all-pairs similarity scoring for protein multiple sequence alignment (MSA). It works on any platform with Vulkan, Metal, or DX12 support and falls back gracefully to CPU when no compatible GPU is available.

---

## Features

| Feature | Description |
|---------|-------------|
| **GPU-accelerated scoring** | Pairwise composition-vector similarity computed on the GPU via a WGSL compute shader |
| **Automatic fallback** | Detects GPU availability at runtime — falls back to CPU seamlessly |
| **Multi-format input** | Reads FASTA, FASTQ, and Stockholm (`.sto`) sequence files |
| **Multi-format output** | Writes aligned FASTA or Clustal W (`.aln`) |
| **Progressive alignment** | Guide-tree-directed Needleman–Wunsch profile alignment |
| **Iterative polishing** | Optional refinement passes improve alignment quality |
| **Dual execution modes** | Run locally with your GPU **or** submit to an HPC cluster via SLURM / PBS |
| **Library + CLI** | Use as a Rust crate (`gpumsa::run`) or as a standalone command-line tool |

---

## Installation

### Prerequisites

- [Rust](https://rustup.rs/) ≥ 1.85 (edition 2024)
- A Vulkan, Metal, or DX12-capable GPU driver (optional — CPU fallback is always available)

### Build from source

```bash
git clone <repository-url>   # replace with your repo URL
cd gpumsa
cargo build --release
```

The binary will be at `target/release/gpumsa`.

### Install as a Cargo binary

```bash
cargo install --path .
```

---

## Usage

### Option 1 — Local GPU

Run directly on your workstation:

```bash
cargo run --release -- \
  --input  examples/toy_proteins.fasta \
  --output results/aligned.fasta \
  --backend auto \
  --refine-passes 1
```

#### CLI Reference

| Flag | Default | Description |
|------|---------|-------------|
| `--input`, `-i` | *(required)* | Input sequence file — FASTA, FASTQ, or Stockholm |
| `--output`, `-o` | *(required)* | Output file — FASTA (default) or Clustal W (`.aln`) |
| `--backend` | `auto` | Similarity backend: `auto`, `gpu`, or `cpu` |
| `--refine-passes` | `1` | Number of iterative polishing passes |
| `--match-score` | `4` | Score for matching residues |
| `--mismatch-score` | `-2` | Penalty for mismatched residues |
| `--gap-penalty` | `-6` | Gap insertion / extension penalty |
| `--check-gpu` | — | Probe the system for a compatible GPU and exit |

#### Supported file formats

**Input** — detected automatically from the file extension:

| Extension | Format |
|-----------|--------|
| `.fa`, `.fasta`, `.faa`, `.fna` | FASTA |
| `.fastq`, `.fq` | FASTQ (quality scores are ignored) |
| `.sto`, `.stk`, `.stockholm` | Stockholm (used by Pfam) |

**Output** — detected automatically from the file extension:

| Extension | Format |
|-----------|--------|
| `.fa`, `.fasta`, `.faa`, `.fna` | FASTA (default) |
| `.aln`, `.clustal` | Clustal W |

#### Examples

```bash
# FASTA input → FASTA output (default)
gpumsa -i examples/toy_proteins.fasta -o results/aligned.fasta

# FASTQ input → FASTA output
gpumsa -i reads.fastq -o results/aligned.fasta

# Stockholm input → Clustal output
gpumsa -i family.sto -o results/aligned.aln

# Force CPU backend, 3 polishing passes
gpumsa -i examples/toy_proteins.fasta -o aligned.fasta \
  --backend cpu --refine-passes 3
```

#### Backend modes

| Mode | Behaviour |
|------|-----------|
| `auto` | Tries GPU first; silently falls back to CPU if no adapter is found |
| `gpu` | Requires a GPU — errors out if no compatible adapter exists |
| `cpu` | CPU-only; useful for deterministic debugging or headless CI |

### Option 2 — HPC Batch (SLURM / PBS)

For large-scale alignments on HPC clusters, use the provided batch scripts.

#### SLURM

```bash
# Default run (uses examples/toy_proteins.fasta)
sbatch scripts/slurm_gpumsa.sh

# Custom input/output via environment variables
INPUT_FASTA=data/pfam_kinases.fasta \
OUTPUT_FASTA=results/kinases_aligned.fasta \
BACKEND=gpu \
REFINE_PASSES=3 \
sbatch scripts/slurm_gpumsa.sh
```

Edit the `#SBATCH` directives at the top of the script to match your cluster:

```bash
#SBATCH --partition=gpu        # your GPU partition name
#SBATCH --gres=gpu:1           # number of GPUs
#SBATCH --cpus-per-task=4      # CPU cores
#SBATCH --mem=16G              # memory
#SBATCH --time=02:00:00        # wall-clock limit
```

#### PBS / Torque

```bash
# Default run
qsub scripts/pbs_gpumsa.sh

# Custom input/output
INPUT_FASTA=data/pfam_kinases.fasta \
OUTPUT_FASTA=results/kinases_aligned.fasta \
qsub scripts/pbs_gpumsa.sh
```

Edit the `#PBS` directives to match your cluster:

```bash
#PBS -q gpu                                       # your GPU queue
#PBS -l select=1:ncpus=4:mem=16gb:ngpus=1         # resources
#PBS -l walltime=02:00:00                         # wall-clock limit
```

#### Environment Variables (both schedulers)

| Variable | Default | Description |
|----------|---------|-------------|
| `INPUT_FASTA` | `examples/toy_proteins.fasta` | Path to input file |
| `OUTPUT_FASTA` | `results/aligned.fasta` | Path to output file |
| `BACKEND` | `auto` | `auto` / `gpu` / `cpu` |
| `REFINE_PASSES` | `1` | Polishing iterations |
| `MATCH_SCORE` | `4` | Match reward |
| `MISMATCH_SCORE` | `-2` | Mismatch penalty |
| `GAP_PENALTY` | `-6` | Gap penalty |

> **Tip:** The batch scripts will automatically build the binary on the compute node if `target/release/gpumsa` does not exist. For large clusters, consider building once on a login node and then submitting jobs.

---

## Architecture

```
gpumsa
├── src/
│   ├── main.rs          # CLI entry point (clap)
│   ├── lib.rs           # Public API re-exports
│   ├── pipeline.rs      # End-to-end orchestration
│   ├── fasta.rs         # Multi-format reader / writer
│   ├── similarity.rs    # Feature extraction + backend dispatch
│   ├── gpu.rs           # wgpu compute shader (WGSL)
│   ├── guide.rs         # Average-linkage guide tree
│   └── aligner.rs       # Needleman–Wunsch progressive alignment
├── examples/
│   ├── align_demo.rs        # Runnable library usage demo
│   ├── toy_proteins.fasta
│   ├── toy_proteins.fastq
│   └── toy_proteins.sto
├── scripts/
│   ├── slurm_gpumsa.sh  # SLURM batch template
│   └── pbs_gpumsa.sh    # PBS batch template
└── Cargo.toml
```

### Pipeline flow

```
Input file (FASTA / FASTQ / Stockholm)
    │
    ▼
┌─────────────────────┐
│  Composition Vectors │   amino-acid frequency features
└─────────┬───────────┘
          │
          ▼
┌─────────────────────┐
│  Pairwise Similarity │   GPU compute shader  ← or CPU fallback
└─────────┬───────────┘
          │
          ▼
┌─────────────────────┐
│  Guide Tree          │   average-linkage clustering
└─────────┬───────────┘
          │
          ▼
┌─────────────────────┐
│  Progressive Align   │   profile–profile Needleman–Wunsch
└─────────┬───────────┘
          │
          ▼
┌─────────────────────┐
│  Polishing           │   iterative re-alignment passes
└─────────┬───────────┘
          │
          ▼
Output file (FASTA / Clustal W)
```

### Module details

| Module | Responsibility |
|--------|----------------|
| `fasta.rs` | Reads FASTA, FASTQ, and Stockholm files (format auto-detected by extension); strips gaps, upper-cases residues. Writes aligned FASTA (80-char line wrapping) or Clustal W. |
| `similarity.rs` | Extracts 21-dimensional amino-acid composition vectors, L2-normalises them, then dispatches to GPU or CPU for all-pairs dot-product scoring. |
| `gpu.rs` | Creates a `wgpu` compute pipeline with a WGSL shader that computes the `n×n` similarity matrix in parallel on the GPU. Uses staging buffers for read-back. |
| `guide.rs` | Builds a binary guide tree from the similarity matrix using greedy average-linkage clustering. |
| `aligner.rs` | Performs progressive profile-profile alignment along the guide tree using Needleman–Wunsch DP, then optionally polishes the result by iteratively re-aligning each sequence against the remaining profile. |
| `pipeline.rs` | Orchestrates the full flow from input to aligned output. Exposes `Config` and `run()` for library consumers. |

---

## Using as a Library

Add to your `Cargo.toml`:

```toml
[dependencies]
gpumsa = { path = "../gpumsa" }  # or publish to crates.io
```

### Public API

Everything you need is re-exported at the crate root — no need to reach into submodules.

| Symbol | Kind | Description |
|--------|------|-------------|
| `run(Config) -> Result<RunSummary>` | fn | Run the full alignment pipeline |
| `Config` | struct | Input path, output path, backend, passes, scoring |
| `RunSummary` | struct | Sequence count, column count, backend used, output path |
| `BackendKind` | enum | `Auto` / `Gpu` / `Cpu` |
| `ScoreScheme` | struct | `match_score`, `mismatch_score`, `gap_penalty` |
| `detect_gpu() -> Option<GpuInfo>` | fn | Probe for a compatible GPU adapter |
| `GpuInfo` | struct | Adapter name, backend, device type, driver |
| `Sequence` | struct | `id: String`, `residues: Vec<u8>` |
| `read_sequences(path)` | fn | Auto-detect format from extension and read |
| `read_fasta(path)` | fn | Read a FASTA file directly |
| `read_fastq(path)` | fn | Read a FASTQ file directly |
| `read_stockholm(path)` | fn | Read a Stockholm file directly |
| `write_alignment(path, &Alignment)` | fn | Write FASTA or Clustal W (by extension) |

### One-liner pipeline

```rust
use gpumsa::{run, Config, BackendKind, ScoreScheme};

let summary = run(Config {
    input: "proteins.fasta".into(),  // also accepts .fastq / .sto
    output: "aligned.fasta".into(),  // use .aln for Clustal W
    backend: BackendKind::Auto,
    refine_passes: 2,
    score_scheme: ScoreScheme { match_score: 4, mismatch_score: -2, gap_penalty: -6 },
})?;

println!("Aligned {} sequences into {} columns ({})",
    summary.sequences, summary.columns, summary.backend);
```

### Working with sequences directly

```rust
use gpumsa::{detect_gpu, read_sequences, write_alignment, Sequence};
use std::path::Path;

// Probe GPU
if let Some(info) = detect_gpu() {
    println!("GPU: {} ({})", info.name, info.backend);
}

// Load sequences — format is auto-detected from the extension
let seqs = read_sequences(Path::new("family.fastq"))?;
println!("{} sequences loaded, first = {}", seqs.len(), seqs[0].id);

// Construct a Sequence manually
let custom = Sequence {
    id: "my_protein".into(),
    residues: b"MKTAYIAKQRQISFVK".to_vec(),
};

// Write an alignment you already have
write_alignment(Path::new("output.aln"), &my_alignment)?;
```

### Running the built-in demo

`examples/align_demo.rs` exercises every public symbol end-to-end:

```bash
cargo run --example align_demo
```

---

## Running Tests

```bash
cargo test
```

41 tests across four modules:

**Parser unit tests (`fasta` module)**

| Test | What it checks |
|------|----------------|
| `parses_multiline_fasta` | Multi-line FASTA records, gap stripping, upper-casing |
| `fasta_strips_gaps_from_prealigned_input` | Pre-aligned FASTA with `-` gaps is stripped before re-alignment |
| `fasta_preserves_ambiguous_amino_acids` | X, B, Z, U, O pass through unchanged |
| `fasta_ignores_blank_lines_between_records` | Blank lines between records are tolerated |
| `fasta_error_on_empty_file` | Empty file → error |
| `fasta_error_on_whitespace_only_file` | Whitespace-only file → error |
| `fasta_error_residues_before_first_header` | Residues before `>` header → error |
| `fasta_error_empty_header_id` | Bare `>` with no ID → error |
| `fasta_error_record_with_no_residues` | Header immediately followed by next header → error |
| `fasta_error_last_record_with_no_residues` | Final record has no residues → error |
| `parses_fastq` | Full 4-line FASTQ records, multi-word description stripped to first word |
| `fastq_strips_gaps_and_uppercases` | FASTQ residue normalisation |
| `fastq_multiword_header_uses_first_word_as_id` | Description after ID is discarded |
| `fastq_error_on_empty_file` | Empty file → error |
| `fastq_error_missing_at_prefix` | Header missing `@` → error |
| `fastq_error_missing_plus_separator` | Missing `+` line → error |
| `fastq_error_truncated_after_sequence` | File ends after sequence line → error |
| `fastq_error_truncated_after_plus` | File ends after `+` → error |
| `parses_stockholm` | Basic Stockholm records, gap/dot removal |
| `stockholm_interleaved_multi_block` | Multi-block interleaved format concatenates residues correctly |
| `stockholm_skips_markup_and_comment_lines` | `#=GF` / `#=GC` lines are ignored |
| `stockholm_error_on_empty_file` | Empty file → error |
| `stockholm_error_only_markup_no_sequences` | Markup-only file → error |
| `stockholm_error_sequences_that_are_only_gaps` | All-gap sequences → error |
| `fasta_output_wraps_at_80_chars` | FASTA output lines are ≤ 80 residue characters |
| `clustal_output_splits_into_60_char_blocks` | Clustal W blocks are ≤ 60 residue columns |
| `write_alignment_errors_on_empty_alignment` | Writing an empty alignment → error |

**Aligner unit tests (`aligner` module)**

| Test | What it checks |
|------|----------------|
| `progressive_alignment_keeps_rows_rectangular` | All rows have the same column count after alignment |
| `polishing_preserves_all_sequences` | Iterative polishing keeps every sequence ID |

**Similarity unit test (`similarity` module)**

| Test | What it checks |
|------|----------------|
| `cpu_similarity_is_symmetric` | CPU similarity matrix is symmetric (`S[i,j] == S[j,i]`) |

**Pipeline integration tests (`pipeline` module)**

| Test | What it checks |
|------|----------------|
| `fasta_input_fasta_output` | End-to-end: FASTA → FASTA |
| `fastq_input_fasta_output` | End-to-end: FASTQ → FASTA |
| `stockholm_input_fasta_output` | End-to-end: Stockholm → FASTA |
| `fasta_input_clustal_output` | End-to-end: FASTA → Clustal W |
| `single_sequence_passes_through_unchanged` | Single sequence aligns without error, output matches input |
| `identical_sequences_align_without_gaps` | Two identical sequences produce zero gaps |
| `unequal_length_sequences_shorter_gets_gaps` | Shorter sequence is padded with `-` gaps |
| `long_sequences_wrap_at_80_in_fasta_output` | 96-residue sequences produce lines ≤ 80 chars |
| `long_sequences_split_into_blocks_in_clustal_output` | 80-residue alignment produces ≥ 2 Clustal blocks |
| `ambiguous_amino_acids_pass_through_pipeline` | X, B, Z, U survive the full pipeline |
| `fastq_and_fasta_inputs_produce_same_alignment` | FASTQ and FASTA files with the same sequences produce byte-identical output |

---

## Limitations

- **O(n²) guide tree** — Memory and runtime scale quadratically with sequence count; not yet tuned for >10k families.
- **CPU-bound progressive stage** — The Needleman–Wunsch DP runs on CPU; only similarity scoring is GPU-accelerated today.
- **Composition-based similarity** — Uses amino-acid frequency vectors rather than substitution matrices or seed-chain methods.
- **No affine gap model** — Uses a simple linear gap penalty; affine (open + extend) is planned.

---

## Roadmap

1. Replace composition similarity with GPU k-mer or seed-diagonal scoring
2. Add a wavefront or striped profile–profile GPU kernel for the progressive stage
3. Stream pair scoring in blocks for >10k sequence scalability
4. Introduce BLOSUM / PAM substitution matrix support
5. Affine gap penalties (open + extend)
6. Benchmark fixtures against BAliBASE, HomFam, and AlphaFold-style family batches

---

## License

This project is provided as-is for research and educational use. See `LICENSE` for details (if present).
