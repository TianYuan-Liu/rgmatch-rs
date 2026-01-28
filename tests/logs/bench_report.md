# Benchmarking Report

## Executive Summary
Benchmarking was conducted to evaluate the performance of `rgmatch-rs` under various thread counts and batch sizes. 
Results indicate that the default batch size (5000) is reasonable but suboptimal; **decreasing** the batch size improves performance and significantly reduces memory usage. 
Additionally, the single-threaded output writer appears to be a throughput bottleneck.

## Batch Size Tuning (8 Threads)
Contrary to the hypothesis that the chunk size might be too small, tests show that **increasing** the chunk size degrades performance and increases memory consumption. **Decreasing** the chunk size improves metrics.

| Batch Size | Time (s) | Peak Memory (MB) | Context Switches (Invol) |
|------------|----------|------------------|--------------------------|
| **500**    | **25.09**| **768**          | 273,438                  |
| 1,000      | 26.62    | 840              | 340,574                  |
| 5,000      | 28.80    | 1,400            | 401,142                  |
| 10,000     | 27.76    | 1,870            | 288,413                  |
| 50,000     | 31.06    | 2,960            | -                        |
| 100,000    | 36.77    | 3,530            | -                        |
| 500,000    | 63.91    | 4,650            | -                        |

**Observation**: Memory usage scales linearly with batch size. Performance degrades with sizes > 10,000. 
**Recommendation**: Reduce default batch size to **500** or **1,000** to lower memory footprint (-45%) and slightly improve speed (+13%).

## Parallelism Efficiency
Comparing single-threaded vs. 8-threaded execution (Batch Size 5000):

| Threads | Time (s) | Speedup | Efficiency |
|---------|----------|---------|------------|
| 1       | 157.42   | 1.0x    | 100%       |
| 8       | 28.80    | 5.46x   | 68%        |

**Analysis**: A 5.46x speedup on 8 cores indicates good but not perfect scalability.

## Identified Bottlenecks

### 1. Output Formatting (Major)
The current implementation performs result formatting (`format_output_line`) in the single **Writer** thread.
- **Volume**: 25.6 million lines generated in ~25-28 seconds.
- **Rate**: ~1 million lines/second.
- **Impact**: The single writer thread is likely CPU-bound formatting strings, causing worker threads to block when the result channel fills up.

**Proposed Fix**: Move `format_output_line` into the Worker threads. Workers should return pre-formatted strings (or byte buffers) to the Writer. The Writer should only be responsible for I/O (`writer.write_all`). This would parallelize the costly string formatting operation.

### 2. Memory Usage with Large Batches
Large batches cause massive memory spikes because 8 concurrent workers each hold a large vector of results. 
- With batch=500,000, memory usage hits 4.6 GB.
- With batch=500, memory usage is only 0.77 GB.

**Fix**: Adopt smaller default batch sizes (e.g., 500-1000).
