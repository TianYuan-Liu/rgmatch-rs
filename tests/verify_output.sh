#!/bin/bash
set -e

PROJECT_DIR="$(pwd)"
DATA_DIR="$PROJECT_DIR/tests/data"
OUTPUT_DIR="$PROJECT_DIR/tests/logs"

mkdir -p "$OUTPUT_DIR"

GTF="$DATA_DIR/full_genome.gtf"
BED="$DATA_DIR/full_peaks.bed"
BINARY="$PROJECT_DIR/target/release/rgmatch"
OUTPUT_1="$OUTPUT_DIR/output_1.txt"
OUTPUT_8="$OUTPUT_DIR/output_8.txt"

if [ ! -f "$GTF" ] || [ ! -f "$BED" ]; then
    echo "Error: Test data not found."
    # Fallback to smaller test files if available or exit
    if [ -f "$DATA_DIR/test.gtf" ] && [ -f "$DATA_DIR/test.bed" ]; then
        echo "Using small test data."
        GTF="$DATA_DIR/test.gtf"
        BED="$DATA_DIR/test.bed"
    else
         echo "No test data found. Exiting."
         exit 1
    fi
fi

echo "Running with 1 thread..."
"$BINARY" -g "$GTF" -b "$BED" -o "$OUTPUT_1" --threads 1

echo "Running with 8 threads..."
"$BINARY" -g "$GTF" -b "$BED" -o "$OUTPUT_8" --threads 8

echo "Comparing outputs..."
if diff "$OUTPUT_1" "$OUTPUT_8"; then
    echo "SUCCESS: Outputs are identical."
else
    echo "FAILURE: Outputs differ."
    exit 1
fi
