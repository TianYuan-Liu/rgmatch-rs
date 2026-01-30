//! BED file parser with gzip support.
//!
//! Parses BED (Browser Extensible Data) files containing genomic regions.

use ahash::AHashMap;
use anyhow::{Context, Result};
use std::fs::File;
use std::io::BufRead;
use std::path::Path;

use crate::parser::util::create_buffered_reader;
use crate::types::Region;

/// Streaming BED file reader for chunked processing.
///
/// This struct provides an iterator-like interface for reading BED files
/// in chunks, enabling memory-efficient processing of large files.
pub struct BedReader {
    reader: Box<dyn BufRead + Send>,
    num_meta_columns: usize,
}

impl BedReader {
    /// Create a new BedReader from a file path (supports .gz).
    pub fn new(path: &Path) -> Result<Self> {
        let file = File::open(path).context("Failed to open BED file")?;
        let reader = create_buffered_reader(file, path);

        Ok(BedReader {
            reader,
            num_meta_columns: 0,
        })
    }

    /// Get the number of metadata columns found so far.
    pub fn num_meta_columns(&self) -> usize {
        self.num_meta_columns
    }

    /// Read the next chunk of regions from the BED file.
    ///
    /// Returns `None` when EOF is reached. The regions are returned in file order,
    /// preserving the original ordering for deterministic output.
    pub fn read_chunk(&mut self, size: usize) -> Result<Option<Vec<Region>>> {
        let mut regions = Vec::with_capacity(size);
        let mut line = String::new();

        while regions.len() < size {
            line.clear();
            let bytes_read = self
                .reader
                .read_line(&mut line)
                .context("Failed to read BED line")?;

            if bytes_read == 0 {
                // EOF reached
                break;
            }

            // Skip empty lines
            let trimmed = line.trim_end();
            if trimmed.is_empty() {
                continue;
            }

            if let Some(region) = self.parse_line(trimmed) {
                regions.push(region);
            }
        }

        if regions.is_empty() {
            Ok(None)
        } else {
            Ok(Some(regions))
        }
    }

    /// Parse a single BED line into a Region.
    fn parse_line(&mut self, line: &str) -> Option<Region> {
        let fields: Vec<&str> = line.split('\t').collect();

        // Need at least 3 columns: chrom, start, end
        if fields.len() < 3 {
            return None;
        }

        let chrom = fields[0].to_string();

        // Try to parse start and end as integers
        // If they fail (e.g., header line), skip this line
        let start: i64 = fields[1].parse().ok()?;
        let end: i64 = fields[2].parse().ok()?;

        // Extract up to 9 additional BED columns as metadata
        let metadata: Vec<String> = fields
            .iter()
            .skip(3)
            .take(9)
            .map(|s| s.to_string())
            .collect();

        // Track the maximum number of metadata columns
        if metadata.len() > self.num_meta_columns {
            self.num_meta_columns = metadata.len();
        }

        Some(Region::new(chrom, start, end, metadata))
    }
}

/// Result of parsing a BED file.
pub struct BedData {
    /// Regions organized by chromosome.
    pub regions_by_chrom: AHashMap<String, Vec<Region>>,
    /// Number of metadata columns found.
    pub num_meta_columns: usize,
}

/// Parse a BED file and return organized region data.
///
/// Supports both plain text and gzip-compressed BED files.
pub fn parse_bed(path: &Path) -> Result<BedData> {
    let file = File::open(path).context("Failed to open BED file")?;
    let reader = create_buffered_reader(file, path);

    parse_bed_reader(reader)
}

/// Parse BED data from a reader.
fn parse_bed_reader<R: BufRead>(reader: R) -> Result<BedData> {
    let mut regions_by_chrom: AHashMap<String, Vec<Region>> = AHashMap::new();
    let mut num_meta_columns = 0;

    for line_result in reader.lines() {
        let line = line_result.context("Failed to read BED line")?;

        // Skip empty lines
        if line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();

        // Need at least 3 columns: chrom, start, end
        if fields.len() < 3 {
            continue;
        }

        let chrom = fields[0].to_string();

        // Try to parse start and end as integers
        // If they fail (e.g., header line), skip this line
        let start: i64 = match fields[1].parse() {
            Ok(v) => v,
            Err(_) => continue, // Skip header lines
        };
        let end: i64 = match fields[2].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };

        // Extract up to 9 additional BED columns as metadata
        let metadata: Vec<String> = fields
            .iter()
            .skip(3)
            .take(9)
            .map(|s| s.to_string())
            .collect();

        // Track the maximum number of metadata columns
        if metadata.len() > num_meta_columns {
            num_meta_columns = metadata.len();
        }

        let region = Region::new(chrom.clone(), start, end, metadata);
        regions_by_chrom.entry(chrom).or_default().push(region);
    }

    Ok(BedData {
        regions_by_chrom,
        num_meta_columns,
    })
}

/// Get standard BED column headers for metadata columns.
pub fn get_bed_headers(num_columns: usize) -> Vec<&'static str> {
    let all_headers = [
        "name",
        "score",
        "strand",
        "thickStart",
        "thickEnd",
        "itemRgb",
        "blockCount",
        "blockSizes",
        "blockStarts",
    ];

    all_headers.iter().take(num_columns).copied().collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    #[test]
    fn test_parse_bed_basic() {
        let bed_content = "chr1\t100\t200\nchrom2\t300\t400\n";

        let reader = BufReader::new(bed_content.as_bytes());
        let result = parse_bed_reader(reader).unwrap();

        assert!(result.regions_by_chrom.contains_key("chr1"));
        assert!(result.regions_by_chrom.contains_key("chrom2"));

        let chr1_regions = &result.regions_by_chrom["chr1"];
        assert_eq!(chr1_regions.len(), 1);
        assert_eq!(chr1_regions[0].start, 100);
        assert_eq!(chr1_regions[0].end, 200);
        assert!(chr1_regions[0].metadata.is_empty());
    }

    #[test]
    fn test_parse_bed_with_metadata() {
        let bed_content = "chr1\t100\t200\tregion1\t500\t+\n";

        let reader = BufReader::new(bed_content.as_bytes());
        let result = parse_bed_reader(reader).unwrap();

        let regions = &result.regions_by_chrom["chr1"];
        assert_eq!(regions[0].metadata.len(), 3);
        assert_eq!(regions[0].metadata[0], "region1");
        assert_eq!(regions[0].metadata[1], "500");
        assert_eq!(regions[0].metadata[2], "+");
        assert_eq!(result.num_meta_columns, 3);
    }

    #[test]
    fn test_parse_bed_skip_header() {
        let bed_content = "chrom\tstart\tend\tname\nchr1\t100\t200\tregion1\n";

        let reader = BufReader::new(bed_content.as_bytes());
        let result = parse_bed_reader(reader).unwrap();

        // Should skip header line (can't parse 'start' as int)
        assert!(result.regions_by_chrom.contains_key("chr1"));
        assert!(!result.regions_by_chrom.contains_key("chrom"));
    }

    #[test]
    fn test_parse_bed_empty_lines() {
        let bed_content = "\nchr1\t100\t200\n\nchr1\t300\t400\n\n";

        let reader = BufReader::new(bed_content.as_bytes());
        let result = parse_bed_reader(reader).unwrap();

        let regions = &result.regions_by_chrom["chr1"];
        assert_eq!(regions.len(), 2);
    }

    #[test]
    fn test_get_bed_headers() {
        assert_eq!(get_bed_headers(0), Vec::<&str>::new());
        assert_eq!(get_bed_headers(3), vec!["name", "score", "strand"]);
        assert_eq!(
            get_bed_headers(9),
            vec![
                "name",
                "score",
                "strand",
                "thickStart",
                "thickEnd",
                "itemRgb",
                "blockCount",
                "blockSizes",
                "blockStarts"
            ]
        );
    }

    #[test]
    fn test_region_methods() {
        let region = Region::new("chr1".to_string(), 100, 200, vec!["test".to_string()]);

        assert_eq!(region.length(), 101);
        assert_eq!(region.midpoint(), 150);
        assert_eq!(region.id(), "chr1_100_200");
    }

    #[test]
    fn test_region_midpoint_integer_division() {
        // Test that midpoint uses integer division
        let region = Region::new("chr1".to_string(), 100, 201, vec![]);
        // (100 + 201) / 2 = 150 (integer division)
        assert_eq!(region.midpoint(), 150);
    }

    #[test]
    fn test_bed_reader_read_chunk() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        // Create a temporary BED file
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "chr1\t100\t200\tregion1").unwrap();
        writeln!(temp_file, "chr1\t300\t400\tregion2").unwrap();
        writeln!(temp_file, "chr2\t500\t600\tregion3").unwrap();
        writeln!(temp_file, "chr2\t700\t800\tregion4").unwrap();
        writeln!(temp_file, "chr3\t900\t1000\tregion5").unwrap();
        temp_file.flush().unwrap();

        let mut reader = BedReader::new(temp_file.path()).unwrap();

        // Read first chunk of 2
        let chunk1 = reader.read_chunk(2).unwrap().unwrap();
        assert_eq!(chunk1.len(), 2);
        assert_eq!(chunk1[0].chrom, "chr1");
        assert_eq!(chunk1[0].start, 100);
        assert_eq!(chunk1[1].chrom, "chr1");
        assert_eq!(chunk1[1].start, 300);

        // Read second chunk of 2
        let chunk2 = reader.read_chunk(2).unwrap().unwrap();
        assert_eq!(chunk2.len(), 2);
        assert_eq!(chunk2[0].chrom, "chr2");
        assert_eq!(chunk2[0].start, 500);

        // Read last chunk (only 1 region left)
        let chunk3 = reader.read_chunk(2).unwrap().unwrap();
        assert_eq!(chunk3.len(), 1);
        assert_eq!(chunk3[0].chrom, "chr3");

        // EOF
        let chunk4 = reader.read_chunk(2).unwrap();
        assert!(chunk4.is_none());

        // Verify metadata columns tracked
        assert_eq!(reader.num_meta_columns(), 1);
    }

    #[test]
    fn test_bed_reader_skips_headers_and_empty_lines() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "chrom\tstart\tend\tname").unwrap(); // header
        writeln!(temp_file).unwrap(); // empty line
        writeln!(temp_file, "chr1\t100\t200\tregion1").unwrap();
        writeln!(temp_file).unwrap(); // empty line
        writeln!(temp_file, "chr1\t300\t400\tregion2").unwrap();
        temp_file.flush().unwrap();

        let mut reader = BedReader::new(temp_file.path()).unwrap();
        let chunk = reader.read_chunk(10).unwrap().unwrap();

        // Should only get 2 valid regions
        assert_eq!(chunk.len(), 2);
        assert_eq!(chunk[0].start, 100);
        assert_eq!(chunk[1].start, 300);
    }
}
