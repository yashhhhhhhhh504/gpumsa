#!/usr/bin/env bash
#===============================================================================
# gpumsa — PBS/Torque batch script for HPC GPU clusters
#
# Usage:
#   qsub scripts/pbs_gpumsa.sh
#
# Customise the PBS directives and the USER CONFIGURATION block below.
#===============================================================================

#PBS -N gpumsa
#PBS -q gpu                        # ← change to your cluster's GPU queue
#PBS -l select=1:ncpus=4:mem=16gb:ngpus=1
#PBS -l walltime=02:00:00
#PBS -o gpumsa_pbs.out
#PBS -e gpumsa_pbs.err
#PBS -j oe

# Move to the directory from which the job was submitted
cd "${PBS_O_WORKDIR:-.}" || exit 1

# ─── USER CONFIGURATION ──────────────────────────────────────────────────────

INPUT_FASTA="${INPUT_FASTA:-examples/toy_proteins.fasta}"
OUTPUT_FASTA="${OUTPUT_FASTA:-results/aligned.fasta}"
BACKEND="${BACKEND:-auto}"          # auto | gpu | cpu
REFINE_PASSES="${REFINE_PASSES:-1}"
MATCH_SCORE="${MATCH_SCORE:-4}"
MISMATCH_SCORE="${MISMATCH_SCORE:-\-2}"
GAP_PENALTY="${GAP_PENALTY:-\-6}"

# ─── ENVIRONMENT SETUP ───────────────────────────────────────────────────────

module purge 2>/dev/null || true
module load rust   2>/dev/null || true
module load cuda   2>/dev/null || true
module load vulkan 2>/dev/null || true

export PATH="$HOME/.cargo/bin:$PATH"

# ─── BUILD ────────────────────────────────────────────────────────────────────
echo "══════════════════════════════════════════════════════════════"
echo " gpumsa — GPU-accelerated Multiple Sequence Aligner"
echo " PBS Job ID   : $PBS_JOBID"
echo " Node         : $(hostname)"
echo " Input        : $INPUT_FASTA"
echo " Output       : $OUTPUT_FASTA"
echo " Backend      : $BACKEND"
echo "══════════════════════════════════════════════════════════════"

if [ ! -f target/release/gpumsa ]; then
    echo "[$(date '+%H:%M:%S')] Building gpumsa in release mode …"
    cargo build --release 2>&1
    if [ $? -ne 0 ]; then
        echo "ERROR: cargo build failed" >&2
        exit 1
    fi
fi

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
