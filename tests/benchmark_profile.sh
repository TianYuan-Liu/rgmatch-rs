#!/bin/bash
# Profile and diagnose bottlenecks in rgmatch-rs
#
# This script runs rgmatch with internal timing instrumentation to determine
# whether the writer thread or worker threads are the bottleneck.
#
# Usage: ./tests/benchmark_profile.sh [threads] [batch_size]
#
# Example: ./tests/benchmark_profile.sh 8 5000

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_DIR"

# Parameters
THREADS="${1:-8}"
BATCH_SIZE="${2:-5000}"

# Test data paths
GTF_FILE="tests/data/full_genome.gtf"
BED_FILE="tests/data/full_peaks.bed"
OUTPUT_DIR="tests/logs"
OUTPUT_FILE="/tmp/profile_output.txt"

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

echo "==================================="
echo "rgmatch-rs Profiling Benchmark"
echo "==================================="
echo ""
echo "Configuration:"
echo "  Threads:    $THREADS"
echo "  Batch size: $BATCH_SIZE"
echo "  GTF file:   $GTF_FILE"
echo "  BED file:   $BED_FILE"
echo ""

# Check if test files exist
if [ ! -f "$GTF_FILE" ]; then
    echo "ERROR: GTF file not found: $GTF_FILE"
    exit 1
fi

if [ ! -f "$BED_FILE" ]; then
    echo "ERROR: BED file not found: $BED_FILE"
    exit 1
fi

# Build with release + frame pointers for samply compatibility
echo "Building with release profile and frame pointers..."
RUSTFLAGS="-C force-frame-pointers=yes" cargo build --release 2>&1

BINARY="./target/release/rgmatch"

if [ ! -f "$BINARY" ]; then
    echo "ERROR: Binary not found: $BINARY"
    exit 1
fi

echo ""
echo "==================================="
echo "Running with internal metrics"
echo "==================================="
echo ""

# Run with timing instrumentation
START_TIME=$(date +%s.%N)

"$BINARY" \
    -g "$GTF_FILE" \
    -b "$BED_FILE" \
    -o "$OUTPUT_FILE" \
    --threads "$THREADS" \
    --batch-size "$BATCH_SIZE" \
    2>&1 | tee "$OUTPUT_DIR/profile_metrics.txt"

END_TIME=$(date +%s.%N)
ELAPSED=$(echo "$END_TIME - $START_TIME" | bc)

echo ""
echo "==================================="
echo "Summary"
echo "==================================="
echo "Total wall-clock time: ${ELAPSED}s"
echo "Output written to: $OUTPUT_FILE"
echo "Metrics saved to: $OUTPUT_DIR/profile_metrics.txt"
echo ""

# Parse key metrics from output
if grep -q "Max pending results:" "$OUTPUT_DIR/profile_metrics.txt"; then
    echo "Key findings from metrics:"
    grep -E "(Max pending|BOTTLENECK|Workers are|Channel)" "$OUTPUT_DIR/profile_metrics.txt" || true
fi

echo ""
echo "==================================="
echo "Next steps for deeper profiling"
echo "==================================="
echo ""
echo "To get a CPU flame graph with samply:"
echo "  cargo install samply"
echo "  samply record $BINARY -g $GTF_FILE -b $BED_FILE -o /tmp/samply_output.txt --threads $THREADS"
echo ""
echo "This will open an interactive flame graph in your browser."
