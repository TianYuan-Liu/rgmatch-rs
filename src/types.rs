//! Core data structures for rgmatch.
//!
//! This module contains the fundamental types used throughout the genomic
//! region-to-gene matching process.

use std::fmt;
use std::str::FromStr;

/// Strand orientation for genomic features.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    Positive,
    Negative,
}

/// Error type for parsing strand from string.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseStrandError;

impl fmt::Display for ParseStrandError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid strand: expected '+' or '-'")
    }
}

impl std::error::Error for ParseStrandError {}

impl FromStr for Strand {
    type Err = ParseStrandError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Strand::Positive),
            "-" => Ok(Strand::Negative),
            _ => Err(ParseStrandError),
        }
    }
}

impl Strand {
    /// Convert strand to string representation.
    pub fn as_str(&self) -> &'static str {
        match self {
            Strand::Positive => "+",
            Strand::Negative => "-",
        }
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

/// Genomic area types for region annotation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Area {
    Tss,
    FirstExon,
    Promoter,
    Tts,
    Intron,
    GeneBody,
    Upstream,
    Downstream,
}

/// Error type for parsing area from string.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseAreaError;

impl fmt::Display for ParseAreaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid area type")
    }
}

impl std::error::Error for ParseAreaError {}

impl FromStr for Area {
    type Err = ParseAreaError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "TSS" => Ok(Area::Tss),
            "1st_EXON" => Ok(Area::FirstExon),
            "PROMOTER" => Ok(Area::Promoter),
            "TTS" => Ok(Area::Tts),
            "INTRON" => Ok(Area::Intron),
            "GENE_BODY" => Ok(Area::GeneBody),
            "UPSTREAM" => Ok(Area::Upstream),
            "DOWNSTREAM" => Ok(Area::Downstream),
            _ => Err(ParseAreaError),
        }
    }
}

impl Area {
    /// Convert area to string representation.
    pub fn as_str(&self) -> &'static str {
        match self {
            Area::Tss => "TSS",
            Area::FirstExon => "1st_EXON",
            Area::Promoter => "PROMOTER",
            Area::Tts => "TTS",
            Area::Intron => "INTRON",
            Area::GeneBody => "GENE_BODY",
            Area::Upstream => "UPSTREAM",
            Area::Downstream => "DOWNSTREAM",
        }
    }
}

impl fmt::Display for Area {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

/// An exon within a transcript.
#[derive(Debug, Clone)]
pub struct Exon {
    pub start: i64,
    pub end: i64,
    /// Exon number within the transcript (set by renumber_exons).
    pub exon_number: Option<String>,
}

impl Exon {
    /// Create a new exon with start and end coordinates.
    pub fn new(start: i64, end: i64) -> Self {
        Exon {
            start,
            end,
            exon_number: None,
        }
    }

    /// Get exon length.
    pub fn length(&self) -> i64 {
        self.end - self.start + 1
    }
}

/// A transcript containing exons.
#[derive(Debug, Clone)]
pub struct Transcript {
    pub transcript_id: String,
    pub exons: Vec<Exon>,
    /// Minimum start coordinate (initialized to i64::MAX).
    pub start: i64,
    /// Maximum end coordinate (initialized to 0).
    pub end: i64,
}

impl Transcript {
    /// Create a new transcript with the given ID.
    pub fn new(transcript_id: String) -> Self {
        Transcript {
            transcript_id,
            exons: Vec::new(),
            start: i64::MAX,
            end: 0,
        }
    }

    /// Add an exon to this transcript.
    pub fn add_exon(&mut self, exon: Exon) {
        self.exons.push(exon);
    }

    /// Set transcript boundaries explicitly.
    pub fn set_length(&mut self, start: i64, end: i64) {
        self.start = start;
        self.end = end;
    }

    /// Calculate transcript boundaries from exon coordinates.
    pub fn calculate_size(&mut self) {
        for exon in &self.exons {
            if exon.start < self.start {
                self.start = exon.start;
            }
            if exon.end > self.end {
                self.end = exon.end;
            }
        }
    }

    /// Renumber exons based on strand orientation.
    ///
    /// Sorts exons by position and assigns exon numbers.
    /// For positive strand: ascending order (1, 2, 3...).
    /// For negative strand: descending order (N, N-1, ...).
    pub fn renumber_exons(&mut self, strand: Strand) {
        // Sort exons by start position
        self.exons.sort_by_key(|e| e.start);

        let n_exons = self.exons.len();

        match strand {
            Strand::Positive => {
                // Positive strand: 1, 2, 3, ...
                for (i, exon) in self.exons.iter_mut().enumerate() {
                    exon.exon_number = Some((i + 1).to_string());
                }
            }
            Strand::Negative => {
                // Negative strand: N, N-1, N-2, ... (reverse numbering)
                for (i, exon) in self.exons.iter_mut().enumerate() {
                    exon.exon_number = Some((n_exons - i).to_string());
                }
            }
        }
    }
}

/// A gene containing transcripts.
#[derive(Debug, Clone)]
pub struct Gene {
    pub gene_id: String,
    pub strand: Strand,
    pub transcripts: Vec<Transcript>,
    /// Minimum start coordinate (initialized to i64::MAX).
    pub start: i64,
    /// Maximum end coordinate (initialized to 0).
    pub end: i64,
}

impl Gene {
    /// Create a new gene with the given ID and strand.
    pub fn new(gene_id: String, strand: Strand) -> Self {
        Gene {
            gene_id,
            strand,
            transcripts: Vec::new(),
            start: i64::MAX,
            end: 0,
        }
    }

    /// Add a transcript to this gene.
    pub fn add_transcript(&mut self, transcript: Transcript) {
        self.transcripts.push(transcript);
    }

    /// Set gene boundaries explicitly.
    pub fn set_length(&mut self, start: i64, end: i64) {
        self.start = start;
        self.end = end;
    }

    /// Calculate gene boundaries from transcript coordinates.
    pub fn calculate_size(&mut self) {
        for transcript in &self.transcripts {
            if transcript.start < self.start {
                self.start = transcript.start;
            }
            if transcript.end > self.end {
                self.end = transcript.end;
            }
        }
    }
}

/// A candidate match between a genomic region and a gene annotation.
#[derive(Debug, Clone)]
pub struct Candidate {
    pub start: i64,
    pub end: i64,
    pub strand: Strand,
    pub exon_number: String,
    pub area: Area,
    pub transcript: String,
    pub gene: String,
    pub distance: i64,
    pub pctg_region: f64,
    pub pctg_area: f64,
    pub tss_distance: i64,
}

impl Candidate {
    /// Create a new candidate.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        start: i64,
        end: i64,
        strand: Strand,
        exon_number: String,
        area: Area,
        transcript: String,
        gene: String,
        distance: i64,
        pctg_region: f64,
        pctg_area: f64,
        tss_distance: i64,
    ) -> Self {
        Candidate {
            start,
            end,
            strand,
            exon_number,
            area,
            transcript,
            gene,
            distance,
            pctg_region,
            pctg_area,
            tss_distance,
        }
    }
}

/// A genomic region from a BED file.
#[derive(Debug, Clone)]
pub struct Region {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub metadata: Vec<String>,
}

impl Region {
    /// Create a new region.
    pub fn new(chrom: String, start: i64, end: i64, metadata: Vec<String>) -> Self {
        Region {
            chrom,
            start,
            end,
            metadata,
        }
    }

    /// Get the region length (end - start + 1).
    pub fn length(&self) -> i64 {
        self.end - self.start + 1
    }

    /// Get the midpoint of the region (integer division).
    pub fn midpoint(&self) -> i64 {
        (self.start + self.end) / 2
    }

    /// Get the region ID (chrom_start_end).
    pub fn id(&self) -> String {
        format!("{}_{}_{}", self.chrom, self.start, self.end)
    }
}

/// Report level for output.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReportLevel {
    Exon,
    Transcript,
    Gene,
}

/// Error type for parsing report level from string.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseReportLevelError;

impl fmt::Display for ParseReportLevelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "invalid report level: expected 'exon', 'transcript', or 'gene'"
        )
    }
}

impl std::error::Error for ParseReportLevelError {}

impl FromStr for ReportLevel {
    type Err = ParseReportLevelError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "exon" => Ok(ReportLevel::Exon),
            "transcript" => Ok(ReportLevel::Transcript),
            "gene" => Ok(ReportLevel::Gene),
            _ => Err(ParseReportLevelError),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strand_parsing() {
        assert_eq!("+".parse::<Strand>(), Ok(Strand::Positive));
        assert_eq!("-".parse::<Strand>(), Ok(Strand::Negative));
        assert!(".".parse::<Strand>().is_err());
    }

    #[test]
    fn test_area_parsing() {
        assert_eq!("TSS".parse::<Area>(), Ok(Area::Tss));
        assert_eq!("1st_EXON".parse::<Area>(), Ok(Area::FirstExon));
        assert_eq!("PROMOTER".parse::<Area>(), Ok(Area::Promoter));
        assert!("INVALID".parse::<Area>().is_err());
    }

    #[test]
    fn test_exon_length() {
        let exon = Exon::new(100, 200);
        assert_eq!(exon.length(), 101);
    }

    #[test]
    fn test_region_midpoint() {
        let region = Region::new("chr1".to_string(), 100, 200, vec![]);
        assert_eq!(region.midpoint(), 150);

        // Test integer division
        let region2 = Region::new("chr1".to_string(), 100, 201, vec![]);
        assert_eq!(region2.midpoint(), 150); // (100 + 201) / 2 = 150 (integer division)
    }

    #[test]
    fn test_transcript_renumber_positive() {
        let mut transcript = Transcript::new("T1".to_string());
        transcript.add_exon(Exon::new(500, 600));
        transcript.add_exon(Exon::new(100, 200));
        transcript.add_exon(Exon::new(300, 400));

        transcript.renumber_exons(Strand::Positive);

        // After sorting and numbering: 100-200 is 1, 300-400 is 2, 500-600 is 3
        assert_eq!(transcript.exons[0].start, 100);
        assert_eq!(transcript.exons[0].exon_number, Some("1".to_string()));
        assert_eq!(transcript.exons[1].start, 300);
        assert_eq!(transcript.exons[1].exon_number, Some("2".to_string()));
        assert_eq!(transcript.exons[2].start, 500);
        assert_eq!(transcript.exons[2].exon_number, Some("3".to_string()));
    }

    #[test]
    fn test_transcript_renumber_negative() {
        let mut transcript = Transcript::new("T1".to_string());
        transcript.add_exon(Exon::new(100, 200));
        transcript.add_exon(Exon::new(300, 400));

        transcript.renumber_exons(Strand::Negative);

        // After sorting: 100-200 first, 300-400 second
        // For negative strand: first (lowest) gets N, last (highest) gets 1
        assert_eq!(transcript.exons[0].start, 100);
        assert_eq!(transcript.exons[0].exon_number, Some("2".to_string()));
        assert_eq!(transcript.exons[1].start, 300);
        assert_eq!(transcript.exons[1].exon_number, Some("1".to_string()));
    }
}
