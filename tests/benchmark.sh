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

# Check if test data exists
if [ ! -f "$GTF" ] || [ ! -f "$BED" ]; then
    echo "Error: Test data not found at $DATA_DIR"
    echo "Expected files: full_genome.gtf, full_peaks.bed"
    exit 1
fi

# ==========================================
# System Info Collection
# ==========================================
echo "Collecting system information..."

OS_VERSION=$(uname -rs)
CPU_CORES=$(sysctl -n hw.ncpu)
TOTAL_MEM_BYTES=$(sysctl -n hw.memsize)
TOTAL_MEM_GB=$(python3 -c "print(round($TOTAL_MEM_BYTES / 1024 / 1024 / 1024, 1))")

# Get chip model (try sysctl first, fallback to system_profiler)
CHIP=$(sysctl -n machdep.cpu.brand_string 2>/dev/null || \
       system_profiler SPHardwareDataType 2>/dev/null | grep "Chip" | awk -F: '{print $2}' | xargs || \
       echo "Unknown")

echo "System: $CHIP, $CPU_CORES cores, $TOTAL_MEM_GB GB RAM"

# ==========================================
# Input File Metrics
# ==========================================
echo "Analyzing input files..."

GTF_SIZE=$(stat -f%z "$GTF")
GTF_SIZE_MB=$(python3 -c "print(round($GTF_SIZE / 1024 / 1024, 1))")
GTF_LINES=$(wc -l < "$GTF" | tr -d ' ')

BED_SIZE=$(stat -f%z "$BED")
BED_SIZE_MB=$(python3 -c "print(round($BED_SIZE / 1024 / 1024, 1))")
BED_LINES=$(wc -l < "$BED" | tr -d ' ')

TOTAL_INPUT_SIZE=$((GTF_SIZE + BED_SIZE))
TOTAL_INPUT_MB=$(python3 -c "print(round($TOTAL_INPUT_SIZE / 1024 / 1024, 1))")

echo "GTF: ${GTF_SIZE_MB} MB (${GTF_LINES} lines)"
echo "BED: ${BED_SIZE_MB} MB (${BED_LINES} lines)"
echo ""

# ==========================================
# Run Benchmarks
# ==========================================
echo "Running benchmark..."
echo "GTF: $GTF"
echo "BED: $BED"
echo ""

# Function to run benchmark and capture metrics
run_benchmark() {
    local threads=$1
    local time_var_name=$2
    local mem_var_name=$3
    local pagefaults_var_name=$4
    local ctx_switches_var_name=$5

    echo "Running with --threads $threads..."

    # Use /usr/bin/time -l to capture memory and other stats
    START=$(python3 -c 'import time; print(time.time())')
    /usr/bin/time -l "$BINARY" -g "$GTF" -b "$BED" -o "$OUTPUT" --threads "$threads" 2> "$TIME_LOG"
    END=$(python3 -c 'import time; print(time.time())')

    local duration=$(python3 -c "print(round($END - $START, 3))")

    # Parse time -l output for memory stats
    # macOS time -l output format varies, so we handle both common formats
    local peak_mem=$(grep -E "maximum resident set size|peak memory footprint" "$TIME_LOG" | head -1 | awk '{print $1}')
    local page_faults=$(grep "page reclaims" "$TIME_LOG" | head -1 | awk '{print $1}' || echo "0")
    local major_faults=$(grep "page faults" "$TIME_LOG" | head -1 | awk '{print $1}' || echo "0")
    local vol_ctx=$(grep "voluntary context switches" "$TIME_LOG" | head -1 | awk '{print $1}' || echo "0")
    local invol_ctx=$(grep "involuntary context switches" "$TIME_LOG" | head -1 | awk '{print $1}' || echo "0")

    # Total context switches
    local total_ctx=$((vol_ctx + invol_ctx))

    # Handle empty values
    [ -z "$peak_mem" ] && peak_mem=0
    [ -z "$page_faults" ] && page_faults=0
    [ -z "$total_ctx" ] && total_ctx=0

    # Export results via global variables (bash doesn't have easy return of multiple values)
    eval "${time_var_name}=$duration"
    eval "${mem_var_name}=$peak_mem"
    eval "${pagefaults_var_name}=$page_faults"
    eval "${ctx_switches_var_name}=$total_ctx"
}

# Run single-threaded benchmark
run_benchmark 1 TIME_1 PEAK_MEM_1 PAGE_FAULTS_1 CTX_SWITCHES_1
OUTPUT_LINES_1=$(wc -l < "$OUTPUT" | tr -d ' ')
echo "Single-threaded time: ${TIME_1}s"

# Run multi-threaded benchmark
run_benchmark 8 TIME_8 PEAK_MEM_8 PAGE_FAULTS_8 CTX_SWITCHES_8
OUTPUT_LINES_8=$(wc -l < "$OUTPUT" | tr -d ' ')
echo "Multi-threaded time: ${TIME_8}s"

# ==========================================
# Calculate Derived Metrics
# ==========================================

# Convert peak memory to MB
PEAK_MEM_MB_1=$(python3 -c "print(round($PEAK_MEM_1 / 1024 / 1024, 1))")
PEAK_MEM_MB_8=$(python3 -c "print(round($PEAK_MEM_8 / 1024 / 1024, 1))")

# Speedup
SPEEDUP=$(python3 -c "print(round($TIME_1 / $TIME_8, 2))")

# Memory overhead percent
MEM_OVERHEAD=$(python3 -c "print(round(($PEAK_MEM_8 - $PEAK_MEM_1) / $PEAK_MEM_1 * 100, 1)) if $PEAK_MEM_1 > 0 else print(0)")

# Throughput metrics
REGIONS_PER_SEC_1=$(python3 -c "print(round($BED_LINES / $TIME_1))")
REGIONS_PER_SEC_8=$(python3 -c "print(round($BED_LINES / $TIME_8))")

THROUGHPUT_MB_S_1=$(python3 -c "print(round($TOTAL_INPUT_MB / $TIME_1, 1))")
THROUGHPUT_MB_S_8=$(python3 -c "print(round($TOTAL_INPUT_MB / $TIME_8, 1))")

LINES_PER_SEC_1=$(python3 -c "print(round($OUTPUT_LINES_1 / $TIME_1))")
LINES_PER_SEC_8=$(python3 -c "print(round($OUTPUT_LINES_8 / $TIME_8))")

# ==========================================
# Generate JSON Output
# ==========================================
cat > "$LOG_DIR/benchmark.json" << EOF
{
  "timestamp": "$(date -u +"%Y-%m-%dT%H:%M:%SZ")",
  "version": "0.1.0",
  "system": {
    "os": "$OS_VERSION",
    "chip": "$CHIP",
    "cpu_cores": $CPU_CORES,
    "total_memory_gb": $TOTAL_MEM_GB
  },
  "input": {
    "gtf_file": "$(basename "$GTF")",
    "gtf_size_mb": $GTF_SIZE_MB,
    "gtf_lines": $GTF_LINES,
    "bed_file": "$(basename "$BED")",
    "bed_size_mb": $BED_SIZE_MB,
    "bed_lines": $BED_LINES
  },
  "runs": [
    {
      "threads": 1,
      "duration_seconds": $TIME_1,
      "peak_memory_mb": $PEAK_MEM_MB_1,
      "page_faults": $PAGE_FAULTS_1,
      "context_switches": $CTX_SWITCHES_1,
      "output_lines": $OUTPUT_LINES_1,
      "regions_per_second": $REGIONS_PER_SEC_1,
      "throughput_mb_per_second": $THROUGHPUT_MB_S_1
    },
    {
      "threads": 8,
      "duration_seconds": $TIME_8,
      "peak_memory_mb": $PEAK_MEM_MB_8,
      "page_faults": $PAGE_FAULTS_8,
      "context_switches": $CTX_SWITCHES_8,
      "output_lines": $OUTPUT_LINES_8,
      "regions_per_second": $REGIONS_PER_SEC_8,
      "throughput_mb_per_second": $THROUGHPUT_MB_S_8
    }
  ],
  "speedup": $SPEEDUP,
  "memory_overhead_percent": $MEM_OVERHEAD
}
EOF

# Cleanup temporary files
rm -f "$OUTPUT" "$TIME_LOG"

# ==========================================
# Console Output
# ==========================================
echo ""
echo "=========================================="
echo "Benchmark Results"
echo "=========================================="
echo "System: $CHIP, $CPU_CORES cores, $TOTAL_MEM_GB GB RAM"
printf "Input:  GTF %.1f GB (%sM lines), BED %.0f MB (%sM lines)\n" \
    "$(python3 -c "print($GTF_SIZE_MB / 1024)")" \
    "$(python3 -c "print(round($GTF_LINES / 1000000, 2))")" \
    "$BED_SIZE_MB" \
    "$(python3 -c "print(round($BED_LINES / 1000000, 2))")"
echo ""
echo "Single-threaded (--threads 1):"
echo "  Time:       ${TIME_1}s"
echo "  Peak Memory: ${PEAK_MEM_MB_1} MB"
printf "  Throughput:  %s regions/s (%.1f MB/s)\n" \
    "$(python3 -c "print('{:,}'.format($REGIONS_PER_SEC_1))")" \
    "$THROUGHPUT_MB_S_1"
echo ""
echo "Multi-threaded (--threads 8):"
echo "  Time:       ${TIME_8}s"
echo "  Peak Memory: ${PEAK_MEM_MB_8} MB (+${MEM_OVERHEAD}%)"
printf "  Throughput:  %s regions/s (%.1f MB/s)\n" \
    "$(python3 -c "print('{:,}'.format($REGIONS_PER_SEC_8))")" \
    "$THROUGHPUT_MB_S_8"
echo ""
echo "Speedup: ${SPEEDUP}x"
echo "=========================================="
echo "Results saved to $LOG_DIR/benchmark.json"
