#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DATA_DIR="$SCRIPT_DIR/data"
LOG_DIR="$SCRIPT_DIR/logs"

mkdir -p "$LOG_DIR"

GTF="$DATA_DIR/full_genome.gtf"
BED="$DATA_DIR/full_peaks.bed"
BINARY="$PROJECT_DIR/target/release/rgmatch"
OUTPUT="$LOG_DIR/bench_output.txt"
TIME_LOG="$LOG_DIR/diag_time.txt"

run_diag() {
    local threads=$1
    local batch_size=$2
    
    echo "================================================="
    echo "Running threads=$threads, batch=$batch_size"
    echo "================================================="
    
    # Use /usr/bin/time -l to capture all stats
    /usr/bin/time -l "$BINARY" -g "$GTF" -b "$BED" -o "$OUTPUT" --threads "$threads" --batch-size "$batch_size" 2> "$TIME_LOG"
    
    cat "$TIME_LOG"
    echo ""
}

echo "Starting Diagnostic Benchmarks..."

run_diag 8 500
run_diag 8 1000
run_diag 8 5000
run_diag 8 10000

# Also run with 1 thread to calculate efficiency
run_diag 1 5000

echo "Done."
