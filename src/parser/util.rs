//! Utility functions for file parsing.

use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Creates a buffered reader that automatically handles gzip-compressed files.
///
/// This function checks if the file path ends with ".gz" and wraps the file
/// in a GzDecoder if so. Otherwise, it returns a plain buffered reader.
pub fn create_buffered_reader(file: File, path: &Path) -> Box<dyn BufRead + Send> {
    if path.to_string_lossy().ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    }
}
