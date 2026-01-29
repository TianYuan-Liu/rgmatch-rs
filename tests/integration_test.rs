use assert_cmd::Command;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[test]
fn test_golden_output_match() -> Result<(), Box<dyn std::error::Error>> {
    let cargo_manifest_dir = env!("CARGO_MANIFEST_DIR");
    let base_dir = Path::new(cargo_manifest_dir);
    // Point to tests/data relative to rgmatch-rs/
    let full_data_dir = base_dir.join("tests").join("data");

    let gtf_path = full_data_dir.join("full_genome.gtf");
    let bed_path = full_data_dir.join("full_peaks.bed");
    let golden_path = full_data_dir.join("full_golden_output_rust.txt");

    // Output file stays in target or tmp, but for now lets put it in tests/data for inspection if it fails
    // actually let's use a temp path or build dir to avoid polluting source tree
    let output_path = base_dir
        .join("tests")
        .join("data")
        .join("integration_test_output.txt");

    // Ensure test data exists
    if !gtf_path.exists() || !bed_path.exists() {
        eprintln!(
            "Test data not found at {:?}. Skipping integration test.",
            full_data_dir
        );
        return Ok(());
    }

    println!("Running integration test with full data...");
    println!("GTF: {:?}", gtf_path);
    println!("BED: {:?}", bed_path);

    // Run the binary
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_rgmatch"));
    cmd.arg("-g")
        .arg(&gtf_path)
        .arg("-b")
        .arg(&bed_path)
        .arg("-o")
        .arg(&output_path)
        .assert()
        .success();

    println!("Binary finished. Comparing outputs...");

    // Open files for streaming comparison
    let output_file = File::open(&output_path)?;
    let golden_file = File::open(&golden_path).map_err(|_| {
        eprintln!(
            "Golden output not found at {:?}. Cannot verify correctness.",
            golden_path
        );
        // If gold missing, we can't verify, but maybe we shouldn't fail if user just wanted to run binary?
        // But this is a test.
        std::io::Error::new(std::io::ErrorKind::NotFound, "Golden file missing")
    })?;

    let output_reader = BufReader::new(output_file);
    let golden_reader = BufReader::new(golden_file);

    // Compare line by line
    for (i, (out_line, gold_line)) in output_reader.lines().zip(golden_reader.lines()).enumerate() {
        let out_line = out_line?;
        let gold_line = gold_line?;

        if out_line != gold_line {
            panic!(
                "Mismatch at line {}: \nExpected: {}\nActual:   {}",
                i + 1,
                gold_line,
                out_line
            );
        }
    }

    // Verify line counts match (if one file is longer)
    // We need to check if either iterator has more lines.
    // zip stops at the shortest.
    // Re-opening or using a slightly different approach would be better but expensive.
    // Let's rely on file size metadata as a quick check first?
    // Or just count lines carefully.

    // Better comparison approach:
    // We already consumed the zipped lines. If sizes differ, we might have missed tail.
    let out_len = std::fs::metadata(&output_path)?.len();
    let gold_len = std::fs::metadata(&golden_path)?.len();

    if out_len != gold_len {
        // Precise check if length differs but content matched so far means one is a prefix of other
        panic!(
            "File sizes differ! Output: {} bytes, Golden: {} bytes. (Content upstream matched)",
            out_len, gold_len
        );
    }

    // Clean up output file on success
    let _ = std::fs::remove_file(output_path);

    Ok(())
}
