#!/usr/bin/env python3
"""Verify Rust-only matches against GTF coordinates."""

import sys

def parse_region(region_str):
    """Parse 'chr_start_end' format."""
    parts = region_str.rsplit('_', 2)
    return {
        'chrom': parts[0],
        'start': int(parts[1]),
        'end': int(parts[2]),
        'midpoint': (int(parts[1]) + int(parts[2])) // 2
    }

def get_gene_coords(gtf_path, gene_id):
    """Extract gene coordinates from GTF."""
    gene_start, gene_end, strand = None, None, None

    with open(gtf_path) as f:
        for line in f:
            if gene_id not in line or line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            if fields[2] == 'gene':
                gene_start = int(fields[3])
                gene_end = int(fields[4])
                strand = fields[6]
                break

    return gene_start, gene_end, strand

def verify_upstream_downstream(region, gene_start, gene_end, reported_distance, config_distance=10000):
    """Verify UPSTREAM/DOWNSTREAM match validity."""
    # Calculate actual distance from region midpoint to nearest gene boundary
    if region['midpoint'] < gene_start:
        actual_distance = gene_start - region['midpoint']
    elif region['midpoint'] > gene_end:
        actual_distance = region['midpoint'] - gene_end
    else:
        actual_distance = 0  # Region overlaps gene

    is_valid = reported_distance <= config_distance
    distance_matches = abs(actual_distance - reported_distance) <= 10  # Allow small rounding

    return {
        'valid': is_valid,
        'reported_distance': reported_distance,
        'calculated_distance': actual_distance,
        'distance_matches': distance_matches,
        'within_threshold': is_valid
    }

def verify_overlap(region, gene_start, gene_end):
    """Verify INTRON/GENE_BODY/FIRST_EXON match validity."""
    # Check if region overlaps gene body
    overlaps = region['start'] < gene_end and region['end'] > gene_start
    return {'valid': overlaps, 'overlaps_gene': overlaps}

def main():
    gtf_path = "tests/data/full_genome.gtf"
    samples_path = "data/analysis/rust_only_sample_main.csv"

    print("Loading samples...")
    with open(samples_path) as f:
        header = f.readline()  # Skip header
        samples = [line.strip().split('\t') for line in f if line.strip()]

    print(f"\nVerifying {len(samples)} samples...\n")

    valid_count = 0
    for i, row in enumerate(samples[:10]):  # Verify first 10
        region_str = row[0]
        gene_id = row[2]
        area = row[5]
        distance = int(row[6]) if row[6] else 0

        region = parse_region(region_str)
        gene_start, gene_end, strand = get_gene_coords(gtf_path, gene_id)

        if gene_start is None:
            print(f"Sample {i+1}: {region_str} -> {gene_id} - GENE NOT FOUND")
            continue

        if area in ('UPSTREAM', 'DOWNSTREAM'):
            result = verify_upstream_downstream(region, gene_start, gene_end, distance)
        else:
            result = verify_overlap(region, gene_start, gene_end)

        verdict = "VALID" if result['valid'] else "INVALID"
        if result['valid']:
            valid_count += 1

        print(f"Sample {i+1}: {region_str}")
        print(f"  Gene: {gene_id} ({gene_start}-{gene_end}, strand={strand})")
        print(f"  Area: {area}, Distance: {distance}")
        print(f"  Result: {result}")
        print(f"  Verdict: {verdict}")
        print()

    print(f"\n=== SUMMARY ===")
    print(f"Valid: {valid_count}/10")

    if valid_count == 10:
        print("\nCONCLUSION: All samples VALID - Rust is finding correct matches that Python missed.")
        print("RECOMMENDATION: Update golden output to match Rust.")
    elif valid_count == 0:
        print("\nCONCLUSION: All samples INVALID - Rust has a lookback bug.")
        print("RECOMMENDATION: Re-implement fine-grained last_index tracking.")
    else:
        print(f"\nCONCLUSION: Mixed results ({valid_count}/10 valid) - deeper investigation needed.")

if __name__ == "__main__":
    main()
