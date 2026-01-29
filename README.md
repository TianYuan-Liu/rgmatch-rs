# rgmatch-rs

A high-performance Rust implementation of the rgmatch tool for genomic interval matching.

## Credits

- **Original Author**: Pedro Furio
- **Current Maintainer (Rust Version)**: Tianyuan Liu

## Overview

rgmatch-rs annotates genomic regions (BED format) with gene features from GTF annotation files. It identifies overlaps with exons, introns, promoters, TSS, TTS, and upstream/downstream regions, outputting detailed annotations for each region-gene association.

## Installation

### From Source

```bash
# Clone the repository
git clone https://github.com/<user>/rgmatch-rs.git
cd rgmatch-rs

# Build in release mode
cargo build --release

# Binary will be at target/release/rgmatch
```

### Requirements

- Rust 1.70 or later
- Cargo (comes with Rust)

## Usage

### Basic Example

```bash
rgmatch -g annotations.gtf.gz -b regions.bed -o output.txt
```

### Command-Line Options

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| `--gtf` | `-g` | Path to GTF annotation file (supports .gz) | Required |
| `--bed` | `-b` | Path to BED file with regions to annotate | Required |
| `--output` | `-o` | Output file path | Required |
| `--report` | `-r` | Report level: `exon`, `transcript`, or `gene` | `exon` |
| `--distance` | `-q` | Maximum distance (kb) for upstream/downstream | `10` |
| `--tss` | `-t` | TSS region size (bp) | `200` |
| `--tts` | `-s` | TTS region size (bp) | `0` |
| `--promoter` | `-p` | Promoter region size (bp) | `1300` |
| `--perc_area` | `-v` | Minimum % of feature covered | `90` |
| `--perc_region` | `-w` | Minimum % of region covered | `50` |
| `--rules` | `-R` | Priority rules (comma-separated) | See below |
| `--gene` | `-G` | GTF attribute for gene ID | `gene_id` |
| `--transcript` | `-T` | GTF attribute for transcript ID | `transcript_id` |
| `--threads` | `-j` | Number of worker threads | `8` |
| `--batch-size` | | Regions per batch | `5000` |

### Report Modes

- **exon** (default): Report at exon level - each exon overlap is reported separately
- **transcript**: Aggregate results by transcript
- **gene**: Aggregate results by gene

### Priority Rules

The `--rules` option controls how annotations are prioritized when a region overlaps multiple features. Default order:

```
TSS,1st_EXON,GENE_BODY,PROMOTER,INTRON,TTS,UPSTREAM,DOWNSTREAM
```

### Examples

```bash
# Basic annotation with default settings
rgmatch -g gencode.v44.gtf.gz -b peaks.bed -o annotated.txt

# Gene-level reporting with custom distance
rgmatch -g gencode.v44.gtf.gz -b peaks.bed -o annotated.txt \
    -r gene -q 50

# Transcript-level with custom TSS/promoter regions
rgmatch -g gencode.v44.gtf.gz -b peaks.bed -o annotated.txt \
    -r transcript -t 500 -p 2000

# Using more threads for large files
rgmatch -g gencode.v44.gtf.gz -b peaks.bed -o annotated.txt -j 16
```

## Output Format

The output is a tab-separated file with columns:

1. Region coordinates and metadata (from BED file)
2. `AREA` - Feature type (TSS, EXON, INTRON, etc.)
3. `GENE` - Gene ID
4. `TRANSCRIPT` - Transcript ID
5. `EXON_NR` - Exon number(s)
6. `STRAND` - Gene strand (+/-)
7. `DISTANCE` - Distance to feature (0 if overlapping)
8. `TSS_DISTANCE` - Distance to transcription start site
9. `PCTG_DHS` - Percentage of region covered
10. `PCTG_AREA` - Percentage of feature covered
11. `EXON_LEN` - Exon length
12. `REGION_LEN` - Region length

## Running Tests

```bash
# Run library unit tests
cargo test --lib

# Run external unit tests
cargo test --test unit_tests

# Run all tests
cargo test
```

Note: Integration tests require large GTF files and are not included in the standard test suite.

## Performance

rgmatch-rs achieves significant speedups over the original Python implementation through:

- Efficient memory layout with sorted gene indices
- Binary search for region-gene matching
- Parallel processing with configurable thread count
- Streaming I/O for large files

## License

MIT License - see [LICENSE](LICENSE) for details.
