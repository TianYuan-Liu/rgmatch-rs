# Unit Test Coverage Task Notes

## Current Status
Unit tests increased from 282 to 362 tests (80 new tests added this iteration).

## Running Tests
```bash
cargo test --test unit_tests  # Unit tests (362 tests)
cargo test --lib              # Library tests (55 tests)
cargo test                    # All tests (~417 total, excluding integration)
```

## Tests Added This Iteration (Iteration 10)

### match_regions_to_genes Extended Tests (8 tests)
- Empty regions/genes handling
- Multiple regions single gene
- Region order preservation
- last_index optimization
- Large gene lengths
- Negative strand genes
- Multiple overlapping genes

### process_candidates_for_output Gene Level Tests (7 tests)
- Single gene single transcript
- Multiple transcripts same area (merging)
- Multiple transcripts different areas
- Multiple genes
- Empty candidates
- Transcript filtering then gene selection
- Merged exon numbers

### Error Type Display Implementation Tests (9 tests)
- ParseStrandError display message
- ParseAreaError display message
- ParseReportLevelError display message
- std::error::Error trait verification
- Clone/Eq trait verification

### Overlap Module Edge Cases (12 tests)
- Exon boundary start/end
- Region in intron
- Region spanning multiple introns
- Single base exon
- Upstream/downstream for positive/negative strand
- TTS enabled scenarios
- Many exons gene
- Distance threshold edge cases

### GTF extract_attribute Extended Tests (8 tests)
- Attribute with equals sign
- Numbers only values
- Empty values
- Values with dots/underscores/hyphens
- Attribute order independence
- No trailing semicolon

### Config Extended Tests (10 tests)
- new() equals default()
- Clone trait
- Debug trait
- max_lookback_distance with tss/tts/promoter/distance dominance
- parse_rules order preservation
- set_distance_kb edge cases

### TssExonInfo/TtsExonInfo Tests (4 tests)
- Field access
- Negative/zero values

### Transcript Level Processing Tests (4 tests)
- Filtering by transcript
- Multiple transcripts
- pctg_region/area filtering

### BED Reader Extended Tests (5 tests)
- Very long lines
- Max metadata columns
- Scientific notation
- Windows line endings
- Mixed whitespace

### check_tss/check_tts Edge Cases (13 tests)
- Single bp regions
- Exact boundary conditions
- Percentage calculations
- Zero distances
- Large zones
- Spanning multiple zones

## Next Steps for Coverage
1. Add integration tests with real BED/GTF sample files
2. Test gzip-compressed file reading (requires test fixtures)
3. Add tests for main.rs CLI argument parsing
4. Add property-based tests for coordinate calculations
5. Consider code coverage analysis with cargo-tarpaulin
