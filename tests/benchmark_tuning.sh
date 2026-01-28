#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DATA_DIR="$SCRIPT_DIR/data"
LOG_DIR="$SCRIPT_DIR/logs"

mkdir -p "$LOG_DIR"

echo "Building release binary..."
cargo build --release --manifest-path="$PROJECT_DIR/Cargo.toml"

GTF="$DATA_DIR/full_genome.gtf"
BED="$DATA_DIR/full_peaks.bed"
BINARY="$PROJECT_DIR/target/release/rgmatch"
OUTPUT="$LOG_DIR/bench_output.txt"
TIME_LOG="$LOG_DIR/time_output.txt"

if [ ! -f "$GTF" ] || [ ! -f "$BED" ]; then
    echo "Error: Test data not found at $DATA_DIR"
    exit 1
fi

run_benchmark_tuning() {
    local threads=$1
    local batch_size=$2
    local label=$3
    
    echo "Running $label (threads=$threads, batch-size=$batch_size)..."
    
    START=$(python3 -c 'import time; print(time.time())')
    /usr/bin/time -l "$BINARY" -g "$GTF" -b "$BED" -o "$OUTPUT" --threads "$threads" --batch-size "$batch_size" 2> "$TIME_LOG"
    END=$(python3 -c 'import time; print(time.time())')
    
    local duration=$(python3 -c "print(round($END - $START, 3))")
    
    local peak_mem=$(grep -E "maximum resident set size|peak memory footprint" "$TIME_LOG" | head -1 | awk '{print $1}')
    [ -z "$peak_mem" ] && peak_mem=0
    
    local peak_mem_mb=$(python3 -c "print(round($peak_mem / 1024 / 1024, 1))")
    
    echo "  Result: Time=${duration}s, Mem=${peak_mem_mb} MB"
}

echo "Starting Tuning Benchmarks..."

# run_benchmark_tuning 8 1000 "Small (1k)"
run_benchmark_tuning 8 5000 "Default (5k)"
run_benchmark_tuning 8 50000 "Large (50k)"
run_benchmark_tuning 8 100000 "Extra Large (100k)"
run_benchmark_tuning 8 500000 "Huge (500k)"

rm -f "$OUTPUT" "$TIME_LOG"
echo "Done."
