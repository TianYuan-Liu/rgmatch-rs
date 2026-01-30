use assert_cmd::Command;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use tempfile::NamedTempFile;

/// Helper function to run a golden test for a specific report level.
///
/// This compares the output of rgmatch at the specified report level
/// against a pre-generated golden output file.
fn run_golden_test(
    report_level: &str,
    golden_filename: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let cargo_manifest_dir = env!("CARGO_MANIFEST_DIR");
    let base_dir = Path::new(cargo_manifest_dir);
    let data_dir = base_dir.join("tests").join("data");

    let gtf_path = data_dir.join("subset_genome.gtf");
    let bed_path = data_dir.join("subset_peaks.bed");
    let golden_path = data_dir.join(golden_filename);

    // Use a temp file for output to avoid polluting source tree
    let output_file = NamedTempFile::new()?;
    let output_path = output_file.path();

    // Run the binary
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_rgmatch"));
    cmd.arg("-g")
        .arg(&gtf_path)
        .arg("-b")
        .arg(&bed_path)
        .arg("-o")
        .arg(output_path)
        .arg("-r")
        .arg(report_level)
        .assert()
        .success();

    // Open files for streaming comparison
    let output_reader = BufReader::new(File::open(output_path)?);
    let golden_reader = BufReader::new(File::open(&golden_path)?);

    // Compare line by line
    let mut line_num = 0;
    for (out_line, gold_line) in output_reader.lines().zip(golden_reader.lines()) {
        line_num += 1;
        let out_line = out_line?;
        let gold_line = gold_line?;

        if out_line != gold_line {
            panic!(
                "Mismatch at line {} (report_level={}): \nExpected: {}\nActual:   {}",
                line_num, report_level, gold_line, out_line
            );
        }
    }

    // Verify file sizes match (to catch if one file has extra lines)
    let out_len = std::fs::metadata(output_path)?.len();
    let gold_len = std::fs::metadata(&golden_path)?.len();

    if out_len != gold_len {
        panic!(
            "File sizes differ for report_level={}: Output: {} bytes, Golden: {} bytes",
            report_level, out_len, gold_len
        );
    }

    Ok(())
}

#[test]
fn test_golden_output_exon() -> Result<(), Box<dyn std::error::Error>> {
    run_golden_test("exon", "subset_golden_output_exon.txt")
}

#[test]
fn test_golden_output_transcript() -> Result<(), Box<dyn std::error::Error>> {
    run_golden_test("transcript", "subset_golden_output_transcript.txt")
}

#[test]
fn test_golden_output_gene() -> Result<(), Box<dyn std::error::Error>> {
    run_golden_test("gene", "subset_golden_output_gene.txt")
}
