# Unit Test Coverage Task Notes

## Current Status
Unit tests increased from 282 to 348 tests (66 new tests added this iteration).

## Tests Added This Iteration (Iteration 10)

### Error Types Display/Error Trait Tests (9 tests)
- ParseStrandError display and Error trait
- ParseAreaError display and Error trait
- ParseReportLevelError display and Error trait
- Clone implementations for all error types

### match_region_to_genes Overlap Cases (12 tests)
- Case 1: Exon before region, region in intron
- Case 2: Exon overlapping partially (left)
- Case 3: Exon completely inside region
- Case 4: Exon overlapping region shifted right
- Case 5: Region completely within exon
- Case 6: Exon totally after region (upstream)
- Negative strand first exon handling
- TTS checking for positive strand
- TSS checking for upstream region
- Multiple transcripts same gene
- Region spanning multiple features

### GtfData Struct Tests (7 tests)
- max_lengths calculation
- Multiple genes max length
- Multiple chromosomes handling
- GtfData clone
- Exon-only GTF (no gene/transcript entries)
- Invalid strand skipping
- Custom gene/transcript tags

### BedReader Edge Cases (9 tests)
- Very long lines (10000 chars)
- Special chromosome names
- Zero coordinates
- Large coordinates (>1 billion)
- Metadata with special characters
- Chunk boundary handling
- Whitespace in metadata
- Negative coordinates

### Config Default Rules Tests (4 tests)
- DEFAULT_RULES order verification
- DEFAULT_RULES length check
- Config rules matches defaults
- Config tags defaults

### match_regions_to_genes Integration (4 tests)
- Empty input handling
- Empty genes handling
- Multiple sorted regions
- max_gene_length usage

### TSS/TTS Boundary Tests (9 tests)
- Exact TSS boundary
- Exact promoter boundary
- Exact TTS boundary
- Zero TSS/promoter values
- Very large values
- Negative strand coordinate flipping
- TTS negative strand no flip

### apply_rules Edge Cases (5 tests)
- All same area different percentages
- Empty rules fallback
- pctg_area then region filtering
- Multiple groups independent
- Order preservation

### select_transcript Edge Cases (4 tests)
- Merge max percentages
- Different areas picks best
- Fallback to first area
- Merged transcript IDs

### Output Formatting Extended (5 tests)
- Varying meta column counts
- All areas in output
- Decimal formatting precision
- Region ID format
- Metadata trimming

## Running Tests
```bash
cargo test --test unit_tests  # Unit tests (348 tests)
cargo test --lib              # Library tests (55 tests)
cargo test                    # All tests (~403 total, excluding integration)
```

## Next Steps for Coverage
1. Add integration tests with real BED/GTF sample files
2. Test gzip-compressed file reading (requires test fixtures)
3. Add tests for main.rs CLI argument parsing
4. Add property-based tests for coordinate calculations
5. Consider code coverage analysis with cargo-tarpaulin
