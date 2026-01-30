# Unit Test Coverage Task Notes

## Current Status
Unit tests increased from 82 to 119 tests (37 new tests added this iteration).

## Tests Added This Iteration (Iteration 4)

### overlap.rs Edge Cases
- `match_region_to_genes_intron_overlap`: Region in intron gap
- `match_region_to_genes_negative_strand`: Negative strand exon numbering
- `match_region_to_genes_upstream_proximity`: Upstream region detection
- `match_region_to_genes_downstream_proximity`: Downstream region detection
- `match_region_to_genes_multiple_genes`: Region overlapping multiple genes
- `match_region_to_genes_gene_body`: Region in gene body (not exons)

### tss.rs Edge Cases
- `tss_exon_info_creation`: Struct creation
- `entirely_in_tss_zone`: 100% TSS overlap
- `entirely_in_promoter_zone`: Promoter-only overlap
- `spans_promoter_and_upstream`: Multi-zone spanning
- `entirely_in_upstream_zone`: Upstream-only
- `neg_strand_entirely_in_tss`: Negative strand TSS
- `neg_strand_upstream`: Negative strand upstream

### tts.rs Edge Cases
- `tts_exon_info_creation`: Struct creation
- `entirely_in_tts_zone`: 100% TTS overlap
- `entirely_in_downstream_zone`: Downstream-only
- `neg_strand_entirely_in_tts`: Negative strand TTS
- `neg_strand_entirely_in_downstream`: Negative strand downstream
- `tts_zero_distance`: Zero TTS distance
- `tts_pctg_calculations`: Percentage verification

### rules.rs Additional Tests
- `empty_grouped_by`: Empty input handling
- `pctg_area_threshold_filter`: pctg_area filtering
- `multiple_groups`: Multiple transcript groups
- `rules_no_matching_area`: No rule matches fallback
- `select_transcript_empty`: Empty input
- `select_transcript_no_rules_match`: Rule fallback
- `select_transcript_multiple_genes`: Multi-gene selection

### output.rs Tests
- Basic output formatting
- Metadata handling (empty, multiple columns)
- Header generation (0, 6, 9 columns)
- Precision verification (2 decimal places)
- Edge cases (zero values, large values, all areas)

## Running Tests
```bash
cargo test --test unit_tests  # Unit tests (119 tests)
cargo test --lib              # Library tests (55 tests)
cargo test                    # All tests (174 total)
```

## Next Steps for Coverage
1. Add tests for `parser/bed.rs` BedReader with gzip files
2. Add tests for `parser/gtf.rs` edge cases (custom ID tags)
3. Add integration tests with real BED/GTF sample files
4. Test Config with various parameter combinations
5. Increase coverage of error paths in parsers
