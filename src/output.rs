//! Output formatting for rgmatch results.
//!
//! This module handles writing formatted output to files with proper
//! column ordering and number formatting.

use anyhow::Result;

use std::io::Write;

use crate::parser::bed::get_bed_headers;
use crate::types::{Candidate, Region};

/// Write the output header.
pub fn write_header<W: Write>(writer: &mut W, num_meta_columns: usize) -> Result<()> {
    let base_header = "Region\tMidpoint\tGene\tTranscript\tExon/Intron\tArea\tDistance\tTSSDistance\tPercRegion\tPercArea";

    if num_meta_columns > 0 {
        let meta_headers = get_bed_headers(num_meta_columns);
        let meta_str = meta_headers.join("\t");
        writeln!(writer, "{}\t{}", base_header, meta_str)?;
    } else {
        writeln!(writer, "{}", base_header)?;
    }

    Ok(())
}

/// Format a single output line for a region-candidate pair.
pub fn format_output_line(region: &Region, candidate: &Candidate) -> String {
    let region_id = region.id();
    let midpoint = region.midpoint();

    // Format percentages with 2 decimal places
    let pctg_region = format!("{:.2}", candidate.pctg_region);
    let pctg_area = format!("{:.2}", candidate.pctg_area);

    // Build base output
    let mut line = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        region_id,
        midpoint,
        candidate.gene,
        candidate.transcript,
        candidate.exon_number,
        candidate.area,
        candidate.distance,
        candidate.tss_distance,
        pctg_region,
        pctg_area
    );

    // Add metadata columns
    if !region.metadata.is_empty() {
        // Join metadata without trailing characters
        let meta_str = region.metadata.join("\t");
        // Remove any trailing newline or whitespace from last element
        let meta_str = meta_str.trim_end();
        line.push('\t');
        line.push_str(meta_str);
    }

    line
}



#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{Area, Strand};

    #[test]
    fn test_format_output_line() {
        let region = Region::new("chr1".to_string(), 100, 200, vec!["name1".to_string()]);
        let candidate = Candidate::new(
            100,
            200,
            Strand::Positive,
            "1".to_string(),
            Area::Tss,
            "T1".to_string(),
            "G1".to_string(),
            50,
            80.123,
            90.456,
            500,
        );

        let line = format_output_line(&region, &candidate);

        assert!(line.contains("chr1_100_200"));
        assert!(line.contains("150")); // midpoint
        assert!(line.contains("G1"));
        assert!(line.contains("T1"));
        assert!(line.contains("TSS"));
        assert!(line.contains("80.12")); // 2 decimal places
        assert!(line.contains("90.46")); // 2 decimal places
        assert!(line.contains("name1"));
    }

    #[test]
    fn test_midpoint_is_integer() {
        // Test that midpoint uses integer division
        let region = Region::new("chr1".to_string(), 100, 201, vec![]);
        let candidate = Candidate::new(
            100,
            200,
            Strand::Positive,
            "1".to_string(),
            Area::Tss,
            "T1".to_string(),
            "G1".to_string(),
            0,
            100.0,
            100.0,
            0,
        );

        let line = format_output_line(&region, &candidate);

        // (100 + 201) / 2 = 150 (integer division)
        assert!(line.contains("\t150\t"));
        assert!(!line.contains("150.5"));
    }

    #[test]
    fn test_negative_pctg_area() {
        let region = Region::new("chr1".to_string(), 100, 200, vec![]);
        let candidate = Candidate::new(
            100,
            200,
            Strand::Positive,
            "1".to_string(),
            Area::Upstream,
            "T1".to_string(),
            "G1".to_string(),
            1000,
            100.0,
            -1.0,
            500,
        );

        let line = format_output_line(&region, &candidate);

        // Should format -1.0 as -1.00
        assert!(line.contains("-1.00"));
    }

    #[test]
    fn test_write_header() {
        let mut output = Vec::new();

        write_header(&mut output, 0).unwrap();
        let header = String::from_utf8(output).unwrap();
        assert!(header.starts_with("Region\tMidpoint\tGene"));
        assert!(!header.contains("name"));

        let mut output = Vec::new();
        write_header(&mut output, 3).unwrap();
        let header = String::from_utf8(output).unwrap();
        assert!(header.contains("name\tscore\tstrand"));
    }
}
