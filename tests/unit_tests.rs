//! Unit tests ported from Python test_unit.py
//!
//! These tests verify the core logic of rgmatch, especially coordinate mirroring
//! and priority rule application.

use rgmatch::config::Config;
use rgmatch::matcher::rules::{apply_rules, select_transcript};
use rgmatch::matcher::tss::{check_tss, TssExonInfo};
use rgmatch::matcher::tts::{check_tts, TtsExonInfo};
use rgmatch::types::{Area, Candidate, Strand, Transcript};

// -------------------------------------------------------------------------
// Helper functions
// -------------------------------------------------------------------------

fn make_candidate(
    area: Area,
    pctg_region: f64,
    pctg_area: f64,
    transcript: &str,
    gene: &str,
    exon_number: &str,
) -> Candidate {
    Candidate::new(
        100,
        200,
        Strand::Positive,
        exon_number.to_string(),
        area,
        transcript.to_string(),
        gene.to_string(),
        0,
        pctg_region,
        pctg_area,
        100,
    )
}

fn default_rules() -> Vec<Area> {
    vec![
        Area::Tss,
        Area::FirstExon,
        Area::Promoter,
        Area::Tts,
        Area::Intron,
        Area::GeneBody,
        Area::Upstream,
        Area::Downstream,
    ]
}

// -------------------------------------------------------------------------
// 1. Coordinate Arithmetic (Fragile Math) Tests
// -------------------------------------------------------------------------

mod test_check_tss {
    use super::*;

    #[test]
    fn test_pos_strand_boundaries() {
        // Exon: [2000, 3000]. TSS @ 2000.
        // TSS zone: [1800, 2000].
        // Promoter zone: [500, 1800) (1300bp long).

        // Case 1: Exactly at TSS boundary (200bp upstream) -> [1800, 1810]
        // 2000 - 1800 = 200. <= 200 TSS.
        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Positive,
            distance: 0,
        };
        let res = check_tss(1800, 1810, &exon, 200.0, 1300.0);
        assert!(
            res.iter().any(|(tag, _, _)| tag == "TSS"),
            "1800 should be TSS: {:?}",
            res
        );

        // Case 2: Just outside TSS boundary -> [1799, 1810]
        // 2000 - 1799 = 201. > 200. Should be PROMOTER.
        let res = check_tss(1799, 1810, &exon, 200.0, 1300.0);
        let tags: Vec<&str> = res.iter().map(|(tag, _, _)| tag.as_str()).collect();
        assert!(tags.contains(&"PROMOTER"));
        assert!(tags.contains(&"TSS"));

        // Case 3: Far upstream (UPSTREAM tag)
        let exon_far = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Positive,
            distance: 1800,
        };
        let res = check_tss(100, 200, &exon_far, 200.0, 1300.0);
        let tags: Vec<&str> = res.iter().map(|(tag, _, _)| tag.as_str()).collect();
        assert!(tags.contains(&"UPSTREAM"));
        assert!(!tags.contains(&"TSS"));
        assert!(!tags.contains(&"PROMOTER"));
    }

    #[test]
    fn test_neg_strand_mirror() {
        // Exon: [2000, 3000]. Strand "-".
        // TSS @ 3000. Upstream > 3000.

        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Negative,
            distance: 0,
        };

        // Case 1: Region [3200, 3210] should be PROMOTER
        let res = check_tss(3200, 3210, &exon, 200.0, 1300.0);
        let tags: Vec<&str> = res.iter().map(|(tag, _, _)| tag.as_str()).collect();
        assert!(tags.contains(&"PROMOTER"));

        // Case 2: TSS Zone Inside [3100, 3150]
        let res = check_tss(3100, 3150, &exon, 200.0, 1300.0);
        assert!(res.iter().any(|(tag, _, _)| tag == "TSS"));
    }

    #[test]
    fn test_integer_division_check_tss() {
        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Positive,
            distance: 0,
        };
        let res = check_tss(1801, 1810, &exon, 200.0, 1300.0);
        assert!(!res.is_empty());
    }

    #[test]
    fn test_zero_length_region_check_tss() {
        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Positive,
            distance: 0,
        };
        let res = check_tss(1900, 1899, &exon, 200.0, 1300.0);
        assert!(res.is_empty());
    }
}

mod test_check_tts {
    use super::*;

    #[test]
    fn test_pos_strand_tts() {
        // Exon [1000, 2000]. TTS @ 2000.
        // Downstream > 2000.
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Positive,
            distance: 0,
        };

        // Case 1: Downstream 100bp [2100, 2150]
        let res = check_tts(2100, 2150, &exon, 200.0);
        assert!(res.iter().any(|(tag, _, _)| tag == "TTS"));
    }

    #[test]
    fn test_neg_strand_tts() {
        // Exon [1000, 2000]. Strand "-".
        // TTS @ 1000 (Start). Downstream < 1000.
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Negative,
            distance: 0,
        };

        // Case 1: Downstream 100bp [850, 900]
        let res = check_tts(850, 900, &exon, 200.0);
        assert!(res.iter().any(|(tag, _, _)| tag == "TTS"));
    }

    #[test]
    fn test_zero_length_region_check_tts() {
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Positive,
            distance: 0,
        };
        let res = check_tts(2100, 2099, &exon, 200.0);
        assert!(res.is_empty());
    }
}

// -------------------------------------------------------------------------
// 2. Filtering Logic Tests
// -------------------------------------------------------------------------

mod test_apply_rules {
    use super::*;
    use ahash::AHashMap;

    #[test]
    fn test_priority_logic() {
        let rules = default_rules();

        let c1 = make_candidate(Area::Intron, 100.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 100.0, 100.0, "T1", "G1", "1"); // Should win
        let c3 = make_candidate(Area::GeneBody, 100.0, 100.0, "T1", "G1", "1");

        let candidates = vec![c1, c2, c3];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("trans1".to_string(), vec![0, 1, 2]);

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);

        assert_eq!(result.len(), 1);
        assert_eq!(result[0].area, Area::Tss);
    }

    #[test]
    fn test_priority_logic_custom_rules() {
        // Change rules order - Intron now higher priority
        let rules = vec![Area::Intron, Area::Tss];

        let c1 = make_candidate(Area::Intron, 100.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 100.0, 100.0, "T1", "G1", "1");

        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("trans1".to_string(), vec![0, 1]);

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);

        assert_eq!(result.len(), 1);
        assert_eq!(result[0].area, Area::Intron);
    }

    #[test]
    fn test_pctg_region_threshold() {
        let rules = vec![Area::Tss, Area::Intron];

        let c1 = make_candidate(Area::Intron, 60.0, 100.0, "T1", "G1", "1"); // Passes
        let c2 = make_candidate(Area::Tss, 40.0, 100.0, "T1", "G1", "1"); // Fails threshold

        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("T1".to_string(), vec![0, 1]);

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);

        assert_eq!(result.len(), 1);
        assert_eq!(result[0].area, Area::Intron);
    }

    #[test]
    fn test_all_below_region_threshold() {
        let rules = vec![Area::Tss, Area::Intron];

        let c1 = make_candidate(Area::Intron, 30.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 40.0, 100.0, "T1", "G1", "1");

        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("T1".to_string(), vec![0, 1]);

        let result = apply_rules(&candidates, &grouped_by, 90.0, 90.0, &rules);

        // Should still pick one based on rules priority
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].area, Area::Tss);
    }

    #[test]
    fn test_max_pctg_region_tiebreaker() {
        let rules = vec![Area::Tss];

        let c1 = make_candidate(Area::Tss, 80.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 90.0, 100.0, "T2", "G1", "1"); // Higher

        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("T1".to_string(), vec![0, 1]);

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);

        assert_eq!(result.len(), 1);
        assert_eq!(result[0].pctg_region, 90.0);
    }

    #[test]
    fn test_same_area_same_pctg_region_tie() {
        let rules = vec![Area::Tss];

        let c1 = make_candidate(Area::Tss, 80.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 80.0, 100.0, "T2", "G1", "1");

        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("T1".to_string(), vec![0, 1]);

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);

        // Both should be reported (tie)
        assert_eq!(result.len(), 2);
    }
}

// -------------------------------------------------------------------------
// 3. Class Logic Tests - Transcript exon numbering
// -------------------------------------------------------------------------

mod test_transcript {
    use super::*;
    use rgmatch::types::Exon;

    #[test]
    fn test_check_exon_numbers_pos() {
        let mut t = Transcript::new("t1".to_string());
        // Add exons in random order
        t.add_exon(Exon::new(500, 600));
        t.add_exon(Exon::new(100, 200));
        t.add_exon(Exon::new(300, 400));

        t.renumber_exons(Strand::Positive);

        // Should be sorted by start: e1, e2, e3
        // Numbering: 1, 2, 3
        assert_eq!(t.exons[0].start, 100);
        assert_eq!(t.exons[0].exon_number, Some("1".to_string()));
        assert_eq!(t.exons[1].start, 300);
        assert_eq!(t.exons[1].exon_number, Some("2".to_string()));
        assert_eq!(t.exons[2].start, 500);
        assert_eq!(t.exons[2].exon_number, Some("3".to_string()));
    }

    #[test]
    fn test_check_exon_numbers_neg() {
        let mut t = Transcript::new("t1".to_string());
        t.add_exon(Exon::new(300, 400));
        t.add_exon(Exon::new(100, 200));

        t.renumber_exons(Strand::Negative);

        // Sorted by start: [100-200, 300-400]
        // For negative strand: first (lowest) gets N, last (highest) gets 1
        assert_eq!(t.exons[0].start, 100);
        assert_eq!(t.exons[0].exon_number, Some("2".to_string()));
        assert_eq!(t.exons[1].start, 300);
        assert_eq!(t.exons[1].exon_number, Some("1".to_string()));
    }
}

// -------------------------------------------------------------------------
// 4. selectTranscript Tests
// -------------------------------------------------------------------------

mod test_select_transcript {
    use super::*;
    use ahash::AHashMap;

    #[test]
    fn test_single_candidate_per_gene() {
        let rules = default_rules();
        let c1 = make_candidate(Area::Tss, 100.0, 100.0, "T1", "G1", "1");
        let candidates = vec![c1];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("G1".to_string(), vec![0]);

        let result = select_transcript(&candidates, &grouped_by, &rules);
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_multiple_candidates_different_areas() {
        let rules = default_rules();
        let c1 = make_candidate(Area::Intron, 100.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 100.0, 100.0, "T2", "G1", "1");
        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("G1".to_string(), vec![0, 1]);

        let result = select_transcript(&candidates, &grouped_by, &rules);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].area, Area::Tss);
    }

    #[test]
    fn test_multiple_candidates_same_area_tie() {
        let rules = default_rules();
        let c1 = make_candidate(Area::Tss, 80.0, 70.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 90.0, 60.0, "T2", "G1", "2");
        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("G1".to_string(), vec![0, 1]);

        let result = select_transcript(&candidates, &grouped_by, &rules);
        assert_eq!(result.len(), 1);
        // Should contain merged transcript info
        assert!(result[0].transcript.contains("T1"));
        assert!(result[0].transcript.contains("T2"));
    }

    #[test]
    fn test_merged_output_pctg_values() {
        let rules = default_rules();
        let c1 = make_candidate(Area::Tss, 80.0, 70.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 90.0, 60.0, "T2", "G1", "3");
        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("G1".to_string(), vec![0, 1]);

        let result = select_transcript(&candidates, &grouped_by, &rules);
        assert_eq!(result[0].pctg_region, 90.0); // max of 80, 90
        assert_eq!(result[0].pctg_area, 70.0); // max of 70, 60
    }

    #[test]
    fn test_merged_exon_numbers() {
        let rules = default_rules();
        let c1 = make_candidate(Area::Tss, 100.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 100.0, 100.0, "T2", "G1", "3");
        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("G1".to_string(), vec![0, 1]);

        let result = select_transcript(&candidates, &grouped_by, &rules);
        assert!(result[0].exon_number.contains("1"));
        assert!(result[0].exon_number.contains("3"));
    }
}

// -------------------------------------------------------------------------
// 5. Config/Rules Tests
// -------------------------------------------------------------------------

mod test_config {
    use super::*;

    #[test]
    fn test_default_config() {
        let config = Config::default();
        assert_eq!(config.rules.len(), 8);
        assert_eq!(config.perc_area, 90.0);
        assert_eq!(config.perc_region, 50.0);
        assert_eq!(config.tss, 200.0);
        assert_eq!(config.tts, 0.0);
        assert_eq!(config.promoter, 1300.0);
        assert_eq!(config.distance, 10000);
    }

    #[test]
    fn test_parse_rules_valid() {
        let mut config = Config::new();
        let result =
            config.parse_rules("DOWNSTREAM,UPSTREAM,GENE_BODY,INTRON,TTS,PROMOTER,1st_EXON,TSS");
        assert!(result);
        assert_eq!(config.rules.len(), 8);
        assert_eq!(config.rules[0], Area::Downstream);
        assert_eq!(config.rules[7], Area::Tss);
    }

    #[test]
    fn test_parse_rules_missing_tags() {
        let mut config = Config::new();
        let result = config.parse_rules("TSS,1st_EXON,PROMOTER");
        assert!(!result);
    }

    #[test]
    fn test_parse_rules_duplicate_tags() {
        let mut config = Config::new();
        let result = config.parse_rules("TSS,TSS,TSS,TSS,TSS,TSS,TSS,TSS");
        assert!(!result);
    }

    #[test]
    fn test_parse_rules_case_sensitive() {
        let mut config = Config::new();
        let result =
            config.parse_rules("tss,1st_exon,promoter,tts,intron,gene_body,upstream,downstream");
        assert!(!result);
    }

    #[test]
    fn test_set_distance_kb() {
        let mut config = Config::new();
        config.set_distance_kb(20);
        assert_eq!(config.distance, 20000);

        config.set_distance_kb(-1);
        assert_eq!(config.distance, 20000); // Should not change
    }
}

// -------------------------------------------------------------------------
// 6. Bug Regression Tests
// -------------------------------------------------------------------------

mod test_bug_regression {
    use super::*;
    use rgmatch::matcher::overlap::match_region_to_genes;
    use rgmatch::types::Exon;
    use rgmatch::{Gene, Region};
    use std::collections::HashSet;

    fn make_test_gene(
        gene_id: &str,
        start: i64,
        end: i64,
        strand: Strand,
        exons: Vec<(i64, i64)>,
    ) -> Gene {
        let mut gene = Gene::new(gene_id.to_string(), strand);
        gene.set_length(start, end);
        let mut transcript = Transcript::new(format!("TRANS_{}", gene_id.replace("GENE", "")));
        for (i, (exon_start, exon_end)) in exons.iter().enumerate() {
            let mut exon = Exon::new(*exon_start, *exon_end);
            exon.exon_number = Some((i + 1).to_string());
            transcript.add_exon(exon);
        }
        transcript.calculate_size();
        transcript.renumber_exons(strand);
        gene.transcripts.push(transcript);
        gene
    }

    /// Bug #1: Test that Case 2 (partial overlap left) doesn't produce duplicate DOWNSTREAM
    #[test]
    fn test_no_duplicate_downstream_case2() {
        // Region [100, 200) partially overlaps exon [51, 150] on the left
        let config = Config::default();
        let region = Region::new("chr1".into(), 100, 200, vec!["region1".into()]);

        // Single-exon gene - triggers Case 2 (partial overlap left on last exon)
        let genes = vec![make_test_gene(
            "GENE001",
            51,
            150,
            Strand::Positive,
            vec![(51, 150)],
        )];

        let last_index = 0;
        let candidates = match_region_to_genes(&region, &genes, &config, last_index);

        // Count DOWNSTREAM candidates for GENE001
        let downstream_count = candidates
            .iter()
            .filter(|c| c.gene == "GENE001" && c.area == Area::Downstream)
            .count();

        assert_eq!(
            downstream_count, 1,
            "GENE001 DOWNSTREAM should appear exactly once, not {}",
            downstream_count
        );
    }

    /// Bug #1: Test that Case 3 (exon inside region) doesn't produce duplicate DOWNSTREAM
    #[test]
    fn test_no_duplicate_downstream_case3() {
        // Region [1000, 1300) completely contains exon [1050, 1200]
        let config = Config::default();
        let region = Region::new("chr1".into(), 1000, 1300, vec!["region2".into()]);

        let genes = vec![make_test_gene(
            "GENE002",
            1050,
            1200,
            Strand::Positive,
            vec![(1050, 1200)],
        )];

        let last_index = 0;
        let candidates = match_region_to_genes(&region, &genes, &config, last_index);

        let downstream_count = candidates
            .iter()
            .filter(|c| c.gene == "GENE002" && c.area == Area::Downstream)
            .count();

        assert_eq!(
            downstream_count, 1,
            "GENE002 DOWNSTREAM should appear exactly once, not {}",
            downstream_count
        );
    }

    /// Bug #2: Test that proximity candidates are preserved when overlapping gene comes later
    #[test]
    fn test_proximity_candidate_preserved() {
        // Region [5000, 5100)
        // GENE003: multi-exon gene ending at 4900 (proximity candidate)
        // GENE004: multi-exon gene with exon overlapping region (overlapping candidate)
        let config = Config::default();
        let region = Region::new("chr1".into(), 5000, 5100, vec!["region3".into()]);

        let genes = vec![
            // GENE003: ends at 4900, 100bp before region - should be DOWNSTREAM proximity
            make_test_gene(
                "GENE003",
                4700,
                4900,
                Strand::Positive,
                vec![(4700, 4750), (4800, 4900)],
            ),
            // GENE004: exon 2 overlaps region at [4950, 5050]
            make_test_gene(
                "GENE004",
                4850,
                5200,
                Strand::Positive,
                vec![(4850, 4900), (4950, 5050)],
            ),
        ];

        let last_index = 0;
        let candidates = match_region_to_genes(&region, &genes, &config, last_index);

        // GENE003 DOWNSTREAM should be preserved (proximity candidate)
        let gene003_downstream = candidates
            .iter()
            .any(|c| c.gene == "GENE003" && c.area == Area::Downstream);

        // GENE004 should also have candidates (overlapping)
        let gene004_present = candidates.iter().any(|c| c.gene == "GENE004");

        assert!(
            gene003_downstream,
            "GENE003 DOWNSTREAM proximity candidate should be preserved"
        );
        assert!(
            gene004_present,
            "GENE004 overlapping candidate should be present"
        );
    }

    /// Combined test: no duplicates across all test regions
    #[test]
    fn test_no_duplicate_lines_overall() {
        let config = Config::default();

        let regions = vec![
            Region::new("chr1".into(), 100, 200, vec!["region1".into()]),
            Region::new("chr1".into(), 1000, 1300, vec!["region2".into()]),
            Region::new("chr1".into(), 5000, 5100, vec!["region3".into()]),
        ];

        let genes = vec![
            make_test_gene("GENE001", 51, 150, Strand::Positive, vec![(51, 150)]),
            make_test_gene("GENE002", 1050, 1200, Strand::Positive, vec![(1050, 1200)]),
            make_test_gene(
                "GENE003",
                4700,
                4900,
                Strand::Positive,
                vec![(4700, 4750), (4800, 4900)],
            ),
            make_test_gene(
                "GENE004",
                4850,
                5200,
                Strand::Positive,
                vec![(4850, 4900), (4950, 5050)],
            ),
        ];

        for region in &regions {
            let last_index = 0;
            let candidates = match_region_to_genes(region, &genes, &config, last_index);

            // Create unique key for each candidate
            let keys: Vec<String> = candidates
                .iter()
                .map(|c| format!("{}_{}_{}", c.gene, c.transcript, c.area))
                .collect();
            let unique_keys: HashSet<_> = keys.iter().collect();

            assert_eq!(
                keys.len(),
                unique_keys.len(),
                "Duplicate candidates found for region {:?}",
                region.id()
            );
        }
    }
}
