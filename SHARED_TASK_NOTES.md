# Unit Test Coverage Task Notes

## Current Status
Unit tests increased from 168 to 226 tests (58 new tests added this iteration).

## Tests Added This Iteration (Iteration 7)

### TSS Extended Tests (10 tests)
- Spanning all zones (TSS, PROMOTER, UPSTREAM)
- Exact boundary tests (positive and negative strand)
- Very small region (1bp)
- Zero TSS/promoter distance handling
- Percentage calculation accuracy

### TTS Extended Tests (8 tests)
- Entirely within TTS zone
- Spanning TTS and DOWNSTREAM
- Negative strand downstream handling
- Exact boundary tests
- Very large TTS zone
- Percentage accuracy

### Rules Extended Tests (9 tests)
- Empty candidates handling
- All fail thresholds fallback
- pctg_area filter behavior
- Multiple groups independence
- Three-candidate merging
- No rules match fallback
- Exact threshold boundary

### Overlap Extended Tests (9 tests)
- find_search_start_index edge cases
- Region completely within exon
- Region spanning multiple exons
- Single exon gene handling
- Beyond distance threshold
- Transcript/gene level processing
- Negative strand first exon

### Output Extended Tests (6 tests)
- Metadata with newlines/whitespace
- Special characters handling
- Exact header output
- Negative coordinates
- Merged transcripts format
- All strands handling

### Parser BED Extended Tests (5 tests)
- Whitespace handling
- Very long lines
- Mixed valid/invalid lines
- Tab-only lines
- Coordinate ordering edge cases

### Parser GTF Extended Tests (7 tests)
- Overlapping genes
- Malformed attributes
- CDS/UTR entries (non-exon)
- Quoted values with spaces
- No exon entries
- Different sources
- Max length multiple chroms

### Config Comprehensive Tests (6 tests)
- Whitespace in rules parsing
- All area combinations
- Extreme values
- Percentage value ranges
- Report level default
- Distance overflow prevention

## Running Tests
```bash
cargo test --test unit_tests  # Unit tests (226 tests)
cargo test --lib              # Library tests (55 tests)
cargo test                    # All tests (~281 total)
```

## Next Steps for Coverage
1. Add integration tests with real BED/GTF sample files
2. Test gzip-compressed file reading (requires test fixtures)
3. Add tests for main.rs CLI argument parsing
4. Add property-based tests for coordinate calculations
5. Consider code coverage analysis with cargo-tarpaulin
