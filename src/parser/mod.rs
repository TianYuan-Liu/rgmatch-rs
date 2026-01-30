//! Parsers for genomic file formats.

pub mod bed;
pub mod gtf;
pub mod util;

pub use bed::{parse_bed, BedReader};
pub use gtf::{parse_gtf, GtfData};
