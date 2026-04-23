#!/usr/bin/env bash
#===============================================================================
# gpumsa — SLURM batch script for HPC GPU clusters
#
# Usage:
#   sbatch scripts/slurm_gpumsa.sh
#
# Customise the SBATCH directives and the USER CONFIGURATION block below.
#===============================================================================

#SBATCH --job-name=gpumsa
#SBATCH --partition=gpu            # ← change to your cluster's GPU partition
#SBATCH --gres=gpu:1               # request 1 GPU
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00            # wall-clock limit
#SBATCH --output=gpumsa_%j.out
#SBATCH --error=gpumsa_%j.err

# ─── USER CONFIGURATION ──────────────────────────────────────────────────────
# Adjust these paths and parameters to match your environment.

INPUT_FASTA="${INPUT_FASTA:-examples/toy_proteins.fasta}"
OUTPUT_FASTA="${OUTPUT_FASTA:-results/aligned.fasta}"
BACKEND="${BACKEND:-auto}"          # auto | gpu | cpu
REFINE_PASSES="${REFINE_PASSES:-1}"
MATCH_SCORE="${MATCH_SCORE:-4}"
MISMATCH_SCORE="${MISMATCH_SCORE:-\-2}"
GAP_PENALTY="${GAP_PENALTY:-\-6}"

# ─── ENVIRONMENT SETUP ───────────────────────────────────────────────────────
# Load modules required on your HPC.  Comment out or modify as needed.
# Common module names are shown; your cluster may differ.

module purge 2>/dev/null || true
module load rust   2>/dev/null || true   # or: module load cargo
module load cuda   2>/dev/null || true   # needed for GPU adapter discovery
module load vulkan 2>/dev/null || true   # wgpu can use Vulkan on Linux

# If Rust is installed via rustup in your home directory, ensure it is on PATH:
export PATH="$HOME/.cargo/bin:$PATH"

# ─── BUILD ────────────────────────────────────────────────────────────────────
echo "══════════════════════════════════════════════════════════════"
echo " gpumsa — GPU-accelerated Multiple Sequence Aligner"
echo " SLURM Job ID : $SLURM_JOB_ID"
echo " Node         : $(hostname)"
echo " GPUs         : ${CUDA_VISIBLE_DEVICES:-none detected}"
echo " Input        : $INPUT_FASTA"
echo " Output       : $OUTPUT_FASTA"
echo " Backend      : $BACKEND"
echo "══════════════════════════════════════════════════════════════"

# Build in release mode (skip if binary already exists)
if [ ! -f target/release/gpumsa ]; then
    echo "[$(date '+%H:%M:%S')] Building gpumsa in release mode …"
    cargo build --release 2>&1
    if [ $? -ne 0 ]; then
        echo "ERROR: cargo build failed" >&2
        exit 1
    fi
fi

# Create output directory if needed
mkdir -p "$(dirname "$OUTPUT_FASTA")"

# ─── RUN ──────────────────────────────────────────────────────────────────────
echo "[$(date '+%H:%M:%S')] Starting alignment …"

./target/release/gpumsa \
    --input  "$INPUT_FASTA" \
    --output "$OUTPUT_FASTA" \
    --backend "$BACKEND" \
    --refine-passes "$REFINE_PASSES" \
    --match-score "$MATCH_SCORE" \
    --mismatch-score "$MISMATCH_SCORE" \
    --gap-penalty "$GAP_PENALTY"

EXIT_CODE=$?

echo "[$(date '+%H:%M:%S')] Finished with exit code $EXIT_CODE"
exit $EXIT_CODE
