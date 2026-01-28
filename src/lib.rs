//! rgmatch - Genomic region-to-gene matching library.
//!
//! This library provides functionality for associating genomic regions (from BED files)
//! with gene annotations (from GTF files) based on positional overlap and proximity.
//!
//! # Features
//!
//! - Parse GTF and BED files (with gzip support)
//! - Match regions to genes considering exon/intron structure
//! - Handle TSS/TTS/promoter regions with strand-aware coordinate transformation
//! - Apply configurable priority rules for tie-breaking
//! - Report at exon, transcript, or gene level
//!
//! # Example
//!
//! ```ignore
//! use rgmatch::config::Config;
//! use rgmatch::parser::{parse_gtf, parse_bed};
//! use rgmatch::matcher::match_regions_to_genes;
//! use rgmatch::output::write_results;
//! use std::path::Path;
//!
//! let config = Config::default();
//! let gtf_data = parse_gtf(Path::new("annotations.gtf"), "gene_id", "transcript_id")?;
//! let bed_data = parse_bed(Path::new("regions.bed"))?;
//!
//! for (chrom, regions) in &bed_data.regions_by_chrom {
//!     if let Some(mut genes) = gtf_data.genes_by_chrom.get(chrom).cloned() {
//!         let results = match_regions_to_genes(&regions, &mut genes, &config);
//!         // Process results...
//!     }
//! }
//! ```

pub mod config;
pub mod matcher;
pub mod output;
pub mod parser;
pub mod types;

pub use config::Config;
pub use parser::{BedReader, GtfData};
pub use types::{Area, Candidate, Gene, Region, ReportLevel, Strand, Transcript};
