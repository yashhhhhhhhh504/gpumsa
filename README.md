# gpumsa

> **GPU-accelerated multiple sequence aligner for large protein families**

`gpumsa` is a Rust library and CLI tool that uses [`wgpu`](https://wgpu.rs) compute shaders to accelerate all-pairs similarity scoring for protein multiple sequence alignment (MSA). It works on any platform with Vulkan, Metal, or DX12 support and falls back gracefully to CPU when no compatible GPU is available.

---

##  Features

| Feature | Description |
|---------|-------------|
| **GPU-accelerated scoring** | Pairwise composition-vector similarity computed on the GPU via a WGSL compute shader |
| **Automatic fallback** | Detects GPU availability at runtime — falls back to CPU seamlessly |
| **Progressive alignment** | Guide-tree-directed Needleman–Wunsch profile alignment |
| **Iterative polishing** | Optional refinement passes improve alignment quality |
| **Dual execution modes** | Run locally with your GPU **or** submit to an HPC cluster via SLURM / PBS |
| **Library + CLI** | Use as a Rust crate (`gpumsa::run`) or as a standalone command-line tool |

---

##  Installation

### Prerequisites

- [Rust](https://rustup.rs/) ≥ 1.85 (edition 2024)
- A Vulkan, Metal, or DX12-capable GPU driver (optional — CPU fallback is always available)

### Build from source

```bash
git clone https://github.com/your-org/gpumsa.git
cd gpumsa
cargo build --release
```

The binary will be at `target/release/gpumsa`.

### Install as a Cargo binary

```bash
cargo install --path .
```

---

##  Usage

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
| `--input`, `-i` | *(required)* | Input FASTA file (protein sequences) |
| `--output`, `-o` | *(required)* | Output aligned FASTA file |
| `--backend` | `auto` | Similarity backend: `auto`, `gpu`, or `cpu` |
| `--refine-passes` | `1` | Number of iterative polishing passes |
| `--match-score` | `4` | Score for matching residues |
| `--mismatch-score` | `-2` | Penalty for mismatched residues |
| `--gap-penalty` | `-6` | Gap insertion / extension penalty |

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
| `INPUT_FASTA` | `examples/toy_proteins.fasta` | Path to input FASTA |
| `OUTPUT_FASTA` | `results/aligned.fasta` | Path to output alignment |
| `BACKEND` | `auto` | `auto` / `gpu` / `cpu` |
| `REFINE_PASSES` | `1` | Polishing iterations |
| `MATCH_SCORE` | `4` | Match reward |
| `MISMATCH_SCORE` | `-2` | Mismatch penalty |
| `GAP_PENALTY` | `-6` | Gap penalty |

> **Tip:** The batch scripts will automatically build the binary on the compute node if `target/release/gpumsa` does not exist. For large clusters, consider building once on a login node and then submitting jobs.

---

##  Architecture

```
gpumsa
├── src/
│   ├── main.rs          # CLI entry point (clap)
│   ├── lib.rs           # Public API re-exports
│   ├── pipeline.rs      # End-to-end orchestration
│   ├── fasta.rs         # FASTA reader / writer
│   ├── similarity.rs    # Feature extraction + backend dispatch
│   ├── gpu.rs           # wgpu compute shader (WGSL)
│   ├── guide.rs         # Average-linkage guide tree
│   └── aligner.rs       # Needleman–Wunsch progressive alignment
├── examples/
│   └── toy_proteins.fasta
├── scripts/
│   ├── slurm_gpumsa.sh  # SLURM batch template
│   └── pbs_gpumsa.sh    # PBS batch template
└── Cargo.toml
```

### Pipeline flow

```
Input FASTA
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
Output aligned FASTA
```

### Module details

| Module | Responsibility |
|--------|----------------|
| `fasta.rs` | Reads multi-record FASTA files, strips gaps, upper-cases residues. Writes aligned FASTA with 80-char line wrapping. |
| `similarity.rs` | Extracts 21-dimensional amino-acid composition vectors, L2-normalises them, then dispatches to GPU or CPU for all-pairs dot-product scoring. |
| `gpu.rs` | Creates a `wgpu` compute pipeline with a WGSL shader that computes the `n×n` similarity matrix in parallel on the GPU. Uses staging buffers for read-back. |
| `guide.rs` | Builds a binary guide tree from the similarity matrix using greedy average-linkage clustering. |
| `aligner.rs` | Performs progressive profile-profile alignment along the guide tree using Needleman–Wunsch DP, then optionally polishes the result by iteratively re-aligning each sequence against the remaining profile. |
| `pipeline.rs` | Orchestrates the full flow from FASTA input to aligned FASTA output. Exposes `Config` and `run()` for library consumers. |

---

##  Using as a Library

Add to your `Cargo.toml`:

```toml
[dependencies]
gpumsa = { path = "../gpumsa" }  # or publish to crates.io
```

```rust
use gpumsa::{run, Config, BackendKind, ScoreScheme};

fn main() -> anyhow::Result<()> {
    let summary = run(Config {
        input: "proteins.fasta".into(),
        output: "aligned.fasta".into(),
        backend: BackendKind::Auto,
        refine_passes: 2,
        score_scheme: ScoreScheme {
            match_score: 4,
            mismatch_score: -2,
            gap_penalty: -6,
        },
    })?;

    println!("Aligned {} sequences into {} columns ({})",
        summary.sequences, summary.columns, summary.backend);
    Ok(())
}
```

---

##  Limitations

- **O(n²) guide tree** — Memory and runtime scale quadratically with sequence count; not yet tuned for >10k families.
- **CPU-bound progressive stage** — The Needleman–Wunsch DP runs on CPU; only similarity scoring is GPU-accelerated today.
- **Composition-based similarity** — Uses amino-acid frequency vectors rather than substitution matrices or seed-chain methods.
- **No affine gap model** — Uses a simple linear gap penalty; affine (open + extend) is planned.

---

##  Roadmap

1. Replace composition similarity with GPU k-mer or seed-diagonal scoring
2. Add a wavefront or striped profile–profile GPU kernel for the progressive stage
3. Stream pair scoring in blocks for >10k sequence scalability
4. Introduce BLOSUM / PAM substitution matrix support
5. Affine gap penalties (open + extend)
6. Benchmark fixtures against BAliBASE, HomFam, and AlphaFold-style family batches

---

##  Running Tests

```bash
cargo test
```

Tests cover FASTA parsing, CPU similarity symmetry, progressive alignment rectangularity, polishing correctness, and the end-to-end CPU pipeline.

---

##  License

This project is provided as-is for research and educational use. See `LICENSE` for details (if present).
