#!/usr/bin/env python3
"""Compare Rust and Python rgmatch outputs, export discrepancies for verification."""

import pandas as pd
import os

# Paths
RUST_OUTPUT = "tests/data/integration_test_output.txt"
GOLDEN_OUTPUT = "tests/data/full_golden_output_rust.txt"
ANALYSIS_DIR = "data/analysis"

os.makedirs(ANALYSIS_DIR, exist_ok=True)

# Column names for rgmatch output
COLUMNS = ["Region", "Midpoint", "Gene", "Transcript", "ExonIntron",
           "Area", "Distance", "TSSDistance", "PercRegion", "PercArea",
           "Col11", "Col12", "Col13"]

def load_output(path):
    """Load rgmatch output file."""
    df = pd.read_csv(path, sep='\t', header=None, names=COLUMNS)
    # Create unique key for matching
    df["_key"] = df["Region"] + "_" + df["Gene"] + "_" + df["Transcript"] + "_" + df["Area"]
    return df

def analyze_pair(df_rs, df_py, label):
    """Analyze differences between Rust and Python outputs."""
    rs_keys = set(df_rs["_key"])
    py_keys = set(df_py["_key"])

    rs_only = rs_keys - py_keys
    py_only = py_keys - rs_keys
    common = rs_keys & py_keys

    print(f"\n=== {label} ===")
    print(f"  Rust total: {len(df_rs)}")
    print(f"  Python total: {len(df_py)}")
    print(f"  Rust-only: {len(rs_only)}")
    print(f"  Python-only: {len(py_only)}")
    print(f"  Common: {len(common)}")

    # Export Rust-only samples for geometric verification
    rust_only_mask = df_rs["_key"].isin(rs_only)
    rust_only_rows = df_rs[rust_only_mask]

    # Sample diverse cases: different Areas, different Genes
    sample_size = min(20, len(rust_only_rows))
    if sample_size > 0:
        # Stratified sample by Area type
        sampled = rust_only_rows.groupby("Area", group_keys=False).apply(
            lambda x: x.head(max(1, sample_size // x.name.__hash__() % 5 + 1))
        ).head(sample_size)

        sample_file = os.path.join(ANALYSIS_DIR, f"rust_only_sample_{label}.csv")
        sampled.to_csv(sample_file, sep="\t", index=False)
        print(f"  Saved {len(sampled)} Rust-only samples to {sample_file}")

    # Analyze which genes contribute most to Rust-only
    if len(rs_only) > 0:
        rust_only_df = df_rs[rust_only_mask]
        gene_counts = rust_only_df["Gene"].value_counts().head(10)
        print(f"\n  Top 10 genes in Rust-only matches:")
        for gene, count in gene_counts.items():
            print(f"    {gene}: {count} lines")

    return {
        "label": label,
        "rust_total": len(df_rs),
        "python_total": len(df_py),
        "rust_only": len(rs_only),
        "python_only": len(py_only),
        "common": len(common)
    }

def main():
    print("Loading outputs...")
    df_rust = load_output(RUST_OUTPUT)
    df_golden = load_output(GOLDEN_OUTPUT)

    results = analyze_pair(df_rust, df_golden, "main")

    # Save summary
    summary_file = os.path.join(ANALYSIS_DIR, "comparison_summary.csv")
    pd.DataFrame([results]).to_csv(summary_file, index=False)
    print(f"\nSaved summary to {summary_file}")

if __name__ == "__main__":
    main()
