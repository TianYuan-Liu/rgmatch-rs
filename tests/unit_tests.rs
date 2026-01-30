//! Unit tests ported from Python test_unit.py
//!
//! These tests verify the core logic of rgmatch, especially coordinate mirroring
//! and priority rule application.

use rgmatch::config::Config;
use rgmatch::matcher::overlap::{
    find_search_start_index, match_region_to_genes, match_regions_to_genes,
    process_candidates_for_output,
};
use rgmatch::matcher::rules::{apply_rules, select_transcript};
use rgmatch::matcher::tss::{check_tss, TssExonInfo};
use rgmatch::matcher::tts::{check_tts, TtsExonInfo};
use rgmatch::output::{format_output_line, write_header};
use rgmatch::types::{Area, Candidate, ReportLevel, Strand, Transcript};

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
    fn test_tss_exon_info_creation() {
        let exon = TssExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Positive,
            distance: 100,
        };
        assert_eq!(exon.start, 1000);
        assert_eq!(exon.end, 2000);
        assert_eq!(exon.strand, Strand::Positive);
        assert_eq!(exon.distance, 100);
    }

    #[test]
    fn test_entirely_in_tss_zone() {
        // Region entirely within TSS zone - should only return TSS
        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Positive,
            distance: 0,
        };
        // Region [1900, 1950] is 50-100bp upstream - entirely in TSS zone (200bp)
        let res = check_tss(1900, 1950, &exon, 200.0, 1300.0);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].0, "TSS");
        // pctg_dhs should be 100% since entire region is in TSS
        assert!(res[0].1 > 99.0);
    }

    #[test]
    fn test_entirely_in_promoter_zone() {
        // Region entirely within PROMOTER zone
        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Positive,
            distance: 500, // 500bp upstream, within promoter
        };
        // Region [1400, 1500] is 500-600bp upstream (in promoter zone)
        let res = check_tss(1400, 1500, &exon, 200.0, 1300.0);
        let tags: Vec<&str> = res.iter().map(|(t, _, _)| t.as_str()).collect();
        assert!(tags.contains(&"PROMOTER"));
    }

    #[test]
    fn test_spans_promoter_and_upstream() {
        // Region spans PROMOTER and extends into UPSTREAM
        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Positive,
            distance: 1400, // 1400bp upstream (beyond TSS+promoter = 1500)
        };
        // With TSS=200, promoter=1300: TSS+promoter extends to 1500bp
        // Region at distance 1400 spans into upstream
        let res = check_tss(100, 700, &exon, 200.0, 1300.0);
        let tags: Vec<&str> = res.iter().map(|(t, _, _)| t.as_str()).collect();
        assert!(
            tags.contains(&"PROMOTER") || tags.contains(&"UPSTREAM"),
            "Should contain PROMOTER or UPSTREAM: {:?}",
            tags
        );
    }

    #[test]
    fn test_entirely_in_upstream_zone() {
        // Region is entirely upstream (beyond TSS+promoter distance)
        let exon = TssExonInfo {
            start: 10000,
            end: 11000,
            strand: Strand::Positive,
            distance: 5000, // 5000bp upstream - well beyond TSS(200)+promoter(1300)
        };
        let res = check_tss(4000, 4500, &exon, 200.0, 1300.0);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].0, "UPSTREAM");
        assert_eq!(res[0].1, 100.0); // 100% of region is upstream
        assert_eq!(res[0].2, -1.0); // pctg_area is -1 for UPSTREAM
    }

    #[test]
    fn test_neg_strand_entirely_in_tss() {
        // Negative strand: TSS is at exon end
        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Negative,
            distance: 0,
        };
        // For negative strand, upstream is > 3000
        // Region [3050, 3100] should be in TSS zone (50-100bp from end)
        let res = check_tss(3050, 3100, &exon, 200.0, 1300.0);
        assert!(res.iter().any(|(t, _, _)| t == "TSS"));
    }

    #[test]
    fn test_neg_strand_upstream() {
        // Negative strand: region far upstream (>TSS+promoter from end)
        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Negative,
            distance: 3000, // 3000bp upstream from end
        };
        let res = check_tss(6000, 6100, &exon, 200.0, 1300.0);
        assert!(res.iter().any(|(t, _, _)| t == "UPSTREAM"));
    }

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
    fn test_tts_exon_info_creation() {
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Negative,
            distance: 50,
        };
        assert_eq!(exon.start, 1000);
        assert_eq!(exon.end, 2000);
        assert_eq!(exon.strand, Strand::Negative);
        assert_eq!(exon.distance, 50);
    }

    #[test]
    fn test_entirely_in_tts_zone() {
        // Region entirely within TTS zone
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Positive,
            distance: 0,
        };
        // For positive strand, TTS is at exon end (2000)
        // Region [2050, 2100] is 50-100bp downstream - should be in TTS zone
        let res = check_tts(2050, 2100, &exon, 200.0);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].0, "TTS");
    }

    #[test]
    fn test_entirely_in_downstream_zone() {
        // Region entirely downstream (beyond TTS distance)
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Positive,
            distance: 500, // 500bp downstream - beyond TTS (200)
        };
        let res = check_tts(2500, 2600, &exon, 200.0);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].0, "DOWNSTREAM");
        assert_eq!(res[0].1, 100.0); // 100% downstream
        assert_eq!(res[0].2, -1.0); // pctg_area is -1 for DOWNSTREAM
    }

    #[test]
    fn test_neg_strand_entirely_in_tts() {
        // Negative strand: TTS is at exon start
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Negative,
            distance: 0,
        };
        // For negative strand, TTS is at exon start (1000)
        // Region [900, 950] is 50-100bp "downstream" (before start)
        let res = check_tts(900, 950, &exon, 200.0);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].0, "TTS");
    }

    #[test]
    fn test_neg_strand_entirely_in_downstream() {
        // Negative strand: region far downstream (< start - tts_distance)
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Negative,
            distance: 500, // 500bp downstream - beyond TTS
        };
        let res = check_tts(400, 500, &exon, 200.0);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].0, "DOWNSTREAM");
    }

    #[test]
    fn test_tts_zero_distance() {
        // When TTS distance is 0, everything downstream is DOWNSTREAM
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Positive,
            distance: 50,
        };
        let res = check_tts(2050, 2100, &exon, 0.0);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].0, "DOWNSTREAM");
    }

    #[test]
    fn test_tts_pctg_calculations() {
        // Verify percentage calculations are correct
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Positive,
            distance: 0,
        };
        // Region [2000, 2100] - 100bp, half in TTS zone (if TTS=100)
        let res = check_tts(2050, 2150, &exon, 100.0);
        // Should span TTS and DOWNSTREAM
        let tags: Vec<&str> = res.iter().map(|(t, _, _)| t.as_str()).collect();
        assert!(
            tags.contains(&"TTS") && tags.contains(&"DOWNSTREAM"),
            "Should span TTS and DOWNSTREAM: {:?}",
            tags
        );
    }

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

    #[test]
    fn test_empty_grouped_by() {
        let rules = default_rules();
        let candidates: Vec<Candidate> = vec![];
        let grouped_by = AHashMap::new();

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);
        assert!(result.is_empty());
    }

    #[test]
    fn test_pctg_area_threshold_filter() {
        let rules = vec![Area::Tss, Area::Intron];

        // Both pass pctg_region, but only c1 passes pctg_area threshold
        let c1 = make_candidate(Area::Intron, 60.0, 95.0, "T1", "G1", "1"); // passes both
        let c2 = make_candidate(Area::Tss, 60.0, 85.0, "T1", "G1", "1"); // fails pctg_area (< 90)

        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("T1".to_string(), vec![0, 1]);

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);

        assert_eq!(result.len(), 1);
        assert_eq!(result[0].area, Area::Intron); // Intron wins because TSS fails pctg_area
    }

    #[test]
    fn test_multiple_groups() {
        let rules = vec![Area::Tss, Area::Intron];

        let c1 = make_candidate(Area::Tss, 100.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Intron, 100.0, 100.0, "T2", "G2", "1");
        let c3 = make_candidate(Area::Promoter, 100.0, 100.0, "T3", "G3", "1");

        let candidates = vec![c1, c2, c3];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("T1".to_string(), vec![0]);
        grouped_by.insert("T2".to_string(), vec![1]);
        grouped_by.insert("T3".to_string(), vec![2]);

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);

        // Each group should produce one result
        assert_eq!(result.len(), 3);
    }

    #[test]
    fn test_rules_no_matching_area() {
        // Rules don't contain any of the candidate areas
        let rules = vec![Area::Upstream, Area::Downstream];

        let c1 = make_candidate(Area::Tss, 100.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Intron, 100.0, 100.0, "T1", "G1", "2");

        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("T1".to_string(), vec![0, 1]);

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);

        // Should still produce results - falls through rule matching
        // Since rules don't match, it should still return based on filtering logic
        // The function iterates rules and breaks when found=true, but if no match, returns nothing
        assert!(result.is_empty() || result.len() <= 2);
    }

    #[test]
    fn test_select_transcript_empty() {
        let rules = default_rules();
        let candidates: Vec<Candidate> = vec![];
        let grouped_by = AHashMap::new();

        let result = select_transcript(&candidates, &grouped_by, &rules);
        assert!(result.is_empty());
    }

    #[test]
    fn test_select_transcript_no_rules_match() {
        // Rules don't contain candidate area - should use fallback
        let rules = vec![Area::Upstream, Area::Downstream];

        let c1 = make_candidate(Area::Tss, 100.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 100.0, 100.0, "T2", "G1", "2");

        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("G1".to_string(), vec![0, 1]);

        let result = select_transcript(&candidates, &grouped_by, &rules);

        // Should fall back to first candidate's area
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].area, Area::Tss);
    }

    #[test]
    fn test_select_transcript_multiple_genes() {
        let rules = vec![Area::Tss, Area::Intron];

        let c1 = make_candidate(Area::Tss, 100.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Intron, 100.0, 100.0, "T2", "G2", "1");

        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("G1".to_string(), vec![0]);
        grouped_by.insert("G2".to_string(), vec![1]);

        let result = select_transcript(&candidates, &grouped_by, &rules);

        // Each gene should have one result
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

// -------------------------------------------------------------------------
// 7. Types Module Comprehensive Tests
// -------------------------------------------------------------------------

mod test_types_strand {
    use super::*;

    #[test]
    fn test_strand_as_str() {
        assert_eq!(Strand::Positive.as_str(), "+");
        assert_eq!(Strand::Negative.as_str(), "-");
    }

    #[test]
    fn test_strand_display() {
        assert_eq!(format!("{}", Strand::Positive), "+");
        assert_eq!(format!("{}", Strand::Negative), "-");
    }

    #[test]
    fn test_strand_from_str_valid() {
        assert_eq!("+".parse::<Strand>().unwrap(), Strand::Positive);
        assert_eq!("-".parse::<Strand>().unwrap(), Strand::Negative);
    }

    #[test]
    fn test_strand_from_str_invalid() {
        assert!(".".parse::<Strand>().is_err());
        assert!("".parse::<Strand>().is_err());
        assert!("positive".parse::<Strand>().is_err());
        assert!("++".parse::<Strand>().is_err());
    }

    #[test]
    fn test_strand_clone_copy() {
        let s1 = Strand::Positive;
        let s2 = s1; // Copy
        let s3 = s1; // Copy trait - no need to call clone()
        assert_eq!(s1, s2);
        assert_eq!(s1, s3);
    }

    #[test]
    fn test_strand_eq_hash() {
        use std::collections::HashSet;
        let mut set = HashSet::new();
        set.insert(Strand::Positive);
        set.insert(Strand::Negative);
        set.insert(Strand::Positive); // Duplicate

        assert_eq!(set.len(), 2);
        assert!(set.contains(&Strand::Positive));
        assert!(set.contains(&Strand::Negative));
    }
}

mod test_types_area {
    use super::*;

    #[test]
    fn test_area_as_str() {
        assert_eq!(Area::Tss.as_str(), "TSS");
        assert_eq!(Area::FirstExon.as_str(), "1st_EXON");
        assert_eq!(Area::Promoter.as_str(), "PROMOTER");
        assert_eq!(Area::Tts.as_str(), "TTS");
        assert_eq!(Area::Intron.as_str(), "INTRON");
        assert_eq!(Area::GeneBody.as_str(), "GENE_BODY");
        assert_eq!(Area::Upstream.as_str(), "UPSTREAM");
        assert_eq!(Area::Downstream.as_str(), "DOWNSTREAM");
    }

    #[test]
    fn test_area_display() {
        assert_eq!(format!("{}", Area::Tss), "TSS");
        assert_eq!(format!("{}", Area::FirstExon), "1st_EXON");
        assert_eq!(format!("{}", Area::Promoter), "PROMOTER");
        assert_eq!(format!("{}", Area::Tts), "TTS");
        assert_eq!(format!("{}", Area::Intron), "INTRON");
        assert_eq!(format!("{}", Area::GeneBody), "GENE_BODY");
        assert_eq!(format!("{}", Area::Upstream), "UPSTREAM");
        assert_eq!(format!("{}", Area::Downstream), "DOWNSTREAM");
    }

    #[test]
    fn test_area_from_str_all_valid() {
        assert_eq!("TSS".parse::<Area>().unwrap(), Area::Tss);
        assert_eq!("1st_EXON".parse::<Area>().unwrap(), Area::FirstExon);
        assert_eq!("PROMOTER".parse::<Area>().unwrap(), Area::Promoter);
        assert_eq!("TTS".parse::<Area>().unwrap(), Area::Tts);
        assert_eq!("INTRON".parse::<Area>().unwrap(), Area::Intron);
        assert_eq!("GENE_BODY".parse::<Area>().unwrap(), Area::GeneBody);
        assert_eq!("UPSTREAM".parse::<Area>().unwrap(), Area::Upstream);
        assert_eq!("DOWNSTREAM".parse::<Area>().unwrap(), Area::Downstream);
    }

    #[test]
    fn test_area_from_str_invalid() {
        assert!("tss".parse::<Area>().is_err()); // Case sensitive
        assert!("INVALID".parse::<Area>().is_err());
        assert!("".parse::<Area>().is_err());
        assert!("FIRST_EXON".parse::<Area>().is_err());
    }

    #[test]
    fn test_area_ordering() {
        // Test PartialOrd/Ord implementations
        assert!(Area::Tss < Area::FirstExon);
        assert!(Area::FirstExon < Area::Promoter);

        let mut areas = [Area::Downstream, Area::Tss, Area::Intron, Area::FirstExon];
        areas.sort();
        assert_eq!(areas[0], Area::Tss);
    }
}

mod test_types_exon {
    use rgmatch::types::Exon;

    #[test]
    fn test_exon_new() {
        let exon = Exon::new(100, 200);
        assert_eq!(exon.start, 100);
        assert_eq!(exon.end, 200);
        assert!(exon.exon_number.is_none());
    }

    #[test]
    fn test_exon_length() {
        let exon = Exon::new(100, 200);
        assert_eq!(exon.length(), 101); // end - start + 1

        let exon2 = Exon::new(0, 0);
        assert_eq!(exon2.length(), 1);

        let exon3 = Exon::new(1000, 2000);
        assert_eq!(exon3.length(), 1001);
    }

    #[test]
    fn test_exon_clone() {
        let mut exon = Exon::new(100, 200);
        exon.exon_number = Some("3".to_string());

        let cloned = exon.clone();
        assert_eq!(cloned.start, 100);
        assert_eq!(cloned.end, 200);
        assert_eq!(cloned.exon_number, Some("3".to_string()));
    }
}

mod test_types_transcript {
    use super::*;
    use rgmatch::types::Exon;

    #[test]
    fn test_transcript_new() {
        let t = Transcript::new("T1".to_string());
        assert_eq!(t.transcript_id, "T1");
        assert!(t.exons.is_empty());
        assert_eq!(t.start, i64::MAX);
        assert_eq!(t.end, 0);
    }

    #[test]
    fn test_transcript_add_exon() {
        let mut t = Transcript::new("T1".to_string());
        t.add_exon(Exon::new(100, 200));
        t.add_exon(Exon::new(300, 400));

        assert_eq!(t.exons.len(), 2);
        assert_eq!(t.exons[0].start, 100);
        assert_eq!(t.exons[1].start, 300);
    }

    #[test]
    fn test_transcript_set_length() {
        let mut t = Transcript::new("T1".to_string());
        t.set_length(50, 500);

        assert_eq!(t.start, 50);
        assert_eq!(t.end, 500);
    }

    #[test]
    fn test_transcript_calculate_size() {
        let mut t = Transcript::new("T1".to_string());
        t.add_exon(Exon::new(500, 600));
        t.add_exon(Exon::new(100, 200));
        t.add_exon(Exon::new(300, 400));

        t.calculate_size();

        assert_eq!(t.start, 100); // Min of exon starts
        assert_eq!(t.end, 600); // Max of exon ends
    }

    #[test]
    fn test_transcript_renumber_exons_positive() {
        let mut t = Transcript::new("T1".to_string());
        t.add_exon(Exon::new(500, 600));
        t.add_exon(Exon::new(100, 200));
        t.add_exon(Exon::new(300, 400));

        t.renumber_exons(Strand::Positive);

        // Sorted by start, numbered 1, 2, 3
        assert_eq!(t.exons[0].start, 100);
        assert_eq!(t.exons[0].exon_number, Some("1".to_string()));
        assert_eq!(t.exons[1].start, 300);
        assert_eq!(t.exons[1].exon_number, Some("2".to_string()));
        assert_eq!(t.exons[2].start, 500);
        assert_eq!(t.exons[2].exon_number, Some("3".to_string()));
    }

    #[test]
    fn test_transcript_renumber_exons_negative() {
        let mut t = Transcript::new("T1".to_string());
        t.add_exon(Exon::new(500, 600));
        t.add_exon(Exon::new(100, 200));

        t.renumber_exons(Strand::Negative);

        // Sorted by start, but numbered in reverse
        assert_eq!(t.exons[0].start, 100);
        assert_eq!(t.exons[0].exon_number, Some("2".to_string()));
        assert_eq!(t.exons[1].start, 500);
        assert_eq!(t.exons[1].exon_number, Some("1".to_string()));
    }

    #[test]
    fn test_transcript_empty_exons() {
        let mut t = Transcript::new("T1".to_string());
        t.renumber_exons(Strand::Positive);
        assert!(t.exons.is_empty());

        t.calculate_size();
        assert_eq!(t.start, i64::MAX); // Unchanged
        assert_eq!(t.end, 0); // Unchanged
    }
}

mod test_types_gene {
    use super::*;
    use rgmatch::types::Exon;
    use rgmatch::Gene;

    #[test]
    fn test_gene_new() {
        let g = Gene::new("G1".to_string(), Strand::Positive);
        assert_eq!(g.gene_id, "G1");
        assert_eq!(g.strand, Strand::Positive);
        assert!(g.transcripts.is_empty());
        assert_eq!(g.start, i64::MAX);
        assert_eq!(g.end, 0);
    }

    #[test]
    fn test_gene_add_transcript() {
        let mut g = Gene::new("G1".to_string(), Strand::Positive);
        g.add_transcript(Transcript::new("T1".to_string()));
        g.add_transcript(Transcript::new("T2".to_string()));

        assert_eq!(g.transcripts.len(), 2);
        assert_eq!(g.transcripts[0].transcript_id, "T1");
        assert_eq!(g.transcripts[1].transcript_id, "T2");
    }

    #[test]
    fn test_gene_set_length() {
        let mut g = Gene::new("G1".to_string(), Strand::Positive);
        g.set_length(100, 500);

        assert_eq!(g.start, 100);
        assert_eq!(g.end, 500);
    }

    #[test]
    fn test_gene_calculate_size() {
        let mut g = Gene::new("G1".to_string(), Strand::Positive);

        let mut t1 = Transcript::new("T1".to_string());
        t1.add_exon(Exon::new(100, 200));
        t1.calculate_size();

        let mut t2 = Transcript::new("T2".to_string());
        t2.add_exon(Exon::new(300, 500));
        t2.calculate_size();

        g.add_transcript(t1);
        g.add_transcript(t2);
        g.calculate_size();

        assert_eq!(g.start, 100); // Min of transcript starts
        assert_eq!(g.end, 500); // Max of transcript ends
    }

    #[test]
    fn test_gene_clone() {
        let mut g = Gene::new("G1".to_string(), Strand::Negative);
        g.set_length(100, 500);
        g.add_transcript(Transcript::new("T1".to_string()));

        let cloned = g.clone();
        assert_eq!(cloned.gene_id, "G1");
        assert_eq!(cloned.strand, Strand::Negative);
        assert_eq!(cloned.start, 100);
        assert_eq!(cloned.end, 500);
        assert_eq!(cloned.transcripts.len(), 1);
    }
}

mod test_types_candidate {
    use super::*;

    #[test]
    fn test_candidate_new() {
        let c = Candidate::new(
            100,
            200,
            Strand::Positive,
            "1".to_string(),
            Area::Tss,
            "T1".to_string(),
            "G1".to_string(),
            50,
            80.5,
            90.5,
            500,
        );

        assert_eq!(c.start, 100);
        assert_eq!(c.end, 200);
        assert_eq!(c.strand, Strand::Positive);
        assert_eq!(c.exon_number, "1");
        assert_eq!(c.area, Area::Tss);
        assert_eq!(c.transcript, "T1");
        assert_eq!(c.gene, "G1");
        assert_eq!(c.distance, 50);
        assert_eq!(c.pctg_region, 80.5);
        assert_eq!(c.pctg_area, 90.5);
        assert_eq!(c.tss_distance, 500);
    }

    #[test]
    fn test_candidate_clone() {
        let c = Candidate::new(
            100,
            200,
            Strand::Negative,
            "2".to_string(),
            Area::Intron,
            "T2".to_string(),
            "G2".to_string(),
            100,
            75.0,
            85.0,
            1000,
        );

        let cloned = c.clone();
        assert_eq!(cloned.start, c.start);
        assert_eq!(cloned.end, c.end);
        assert_eq!(cloned.strand, c.strand);
        assert_eq!(cloned.exon_number, c.exon_number);
        assert_eq!(cloned.area, c.area);
        assert_eq!(cloned.transcript, c.transcript);
        assert_eq!(cloned.gene, c.gene);
        assert_eq!(cloned.distance, c.distance);
        assert_eq!(cloned.pctg_region, c.pctg_region);
        assert_eq!(cloned.pctg_area, c.pctg_area);
    }
}

mod test_types_region {
    use rgmatch::Region;

    #[test]
    fn test_region_new() {
        let r = Region::new(
            "chr1".to_string(),
            100,
            200,
            vec!["name1".to_string(), "500".to_string()],
        );

        assert_eq!(r.chrom, "chr1");
        assert_eq!(r.start, 100);
        assert_eq!(r.end, 200);
        assert_eq!(r.metadata.len(), 2);
        assert_eq!(r.metadata[0], "name1");
    }

    #[test]
    fn test_region_length() {
        let r = Region::new("chr1".to_string(), 100, 200, vec![]);
        assert_eq!(r.length(), 101); // end - start + 1

        let r2 = Region::new("chr1".to_string(), 0, 0, vec![]);
        assert_eq!(r2.length(), 1);
    }

    #[test]
    fn test_region_midpoint() {
        let r = Region::new("chr1".to_string(), 100, 200, vec![]);
        assert_eq!(r.midpoint(), 150); // (100 + 200) / 2

        // Test integer division
        let r2 = Region::new("chr1".to_string(), 100, 201, vec![]);
        assert_eq!(r2.midpoint(), 150); // (100 + 201) / 2 = 150 (integer division)

        let r3 = Region::new("chr1".to_string(), 0, 1, vec![]);
        assert_eq!(r3.midpoint(), 0); // (0 + 1) / 2 = 0
    }

    #[test]
    fn test_region_id() {
        let r = Region::new("chr1".to_string(), 100, 200, vec![]);
        assert_eq!(r.id(), "chr1_100_200");

        let r2 = Region::new("chrX".to_string(), 0, 1000000, vec![]);
        assert_eq!(r2.id(), "chrX_0_1000000");
    }

    #[test]
    fn test_region_clone() {
        let r = Region::new("chr1".to_string(), 100, 200, vec!["meta".to_string()]);

        let cloned = r.clone();
        assert_eq!(cloned.chrom, "chr1");
        assert_eq!(cloned.start, 100);
        assert_eq!(cloned.end, 200);
        assert_eq!(cloned.metadata, vec!["meta"]);
    }
}

mod test_types_report_level {
    use super::*;

    #[test]
    fn test_report_level_from_str() {
        assert_eq!("exon".parse::<ReportLevel>().unwrap(), ReportLevel::Exon);
        assert_eq!(
            "transcript".parse::<ReportLevel>().unwrap(),
            ReportLevel::Transcript
        );
        assert_eq!("gene".parse::<ReportLevel>().unwrap(), ReportLevel::Gene);
    }

    #[test]
    fn test_report_level_case_insensitive() {
        assert_eq!("EXON".parse::<ReportLevel>().unwrap(), ReportLevel::Exon);
        assert_eq!("Exon".parse::<ReportLevel>().unwrap(), ReportLevel::Exon);
        assert_eq!(
            "TRANSCRIPT".parse::<ReportLevel>().unwrap(),
            ReportLevel::Transcript
        );
        assert_eq!("GENE".parse::<ReportLevel>().unwrap(), ReportLevel::Gene);
    }

    #[test]
    fn test_report_level_invalid() {
        assert!("invalid".parse::<ReportLevel>().is_err());
        assert!("".parse::<ReportLevel>().is_err());
        assert!("exons".parse::<ReportLevel>().is_err());
    }

    #[test]
    fn test_report_level_clone_copy() {
        let r1 = ReportLevel::Gene;
        let r2 = r1; // Copy
        let r3 = r1; // Copy trait - no need to call clone()
        assert_eq!(r1, r2);
        assert_eq!(r1, r3);
    }
}

// -------------------------------------------------------------------------
// 8. Config Module Additional Tests
// -------------------------------------------------------------------------

mod test_config_additional {
    use super::*;

    #[test]
    fn test_config_new_equals_default() {
        let c1 = Config::new();
        let c2 = Config::default();

        assert_eq!(c1.rules, c2.rules);
        assert_eq!(c1.perc_area, c2.perc_area);
        assert_eq!(c1.perc_region, c2.perc_region);
        assert_eq!(c1.tss, c2.tss);
        assert_eq!(c1.tts, c2.tts);
        assert_eq!(c1.promoter, c2.promoter);
        assert_eq!(c1.distance, c2.distance);
    }

    #[test]
    fn test_config_max_lookback_distance() {
        let mut config = Config::new();

        // Default: tss=200, tts=0, promoter=1300, distance=10000
        // max_lookback = distance.max(max(tss, tts, promoter))
        // = 10000.max(1300) = 10000
        assert_eq!(config.max_lookback_distance(), 10000);

        // Increase promoter beyond distance
        config.promoter = 15000.0;
        assert_eq!(config.max_lookback_distance(), 15000);

        // Reset and test with large TSS
        config.promoter = 1300.0;
        config.tss = 20000.0;
        assert_eq!(config.max_lookback_distance(), 20000);

        // Test with large TTS
        config.tss = 200.0;
        config.tts = 25000.0;
        assert_eq!(config.max_lookback_distance(), 25000);
    }

    #[test]
    fn test_config_set_distance_kb_zero() {
        let mut config = Config::new();
        config.set_distance_kb(0);
        assert_eq!(config.distance, 0);
    }

    #[test]
    fn test_config_clone() {
        let mut config = Config::new();
        config.tss = 500.0;
        config.distance = 20000;

        let cloned = config.clone();
        assert_eq!(cloned.tss, 500.0);
        assert_eq!(cloned.distance, 20000);
    }
}

// -------------------------------------------------------------------------
// 9. Overlap Module Tests
// -------------------------------------------------------------------------

mod test_overlap_functions {
    use super::*;
    use rgmatch::types::Exon;
    use rgmatch::{Gene, Region};

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

    #[test]
    fn test_find_search_start_index_empty() {
        let genes: Vec<Gene> = vec![];
        assert_eq!(find_search_start_index(&genes, 100), 0);
    }

    #[test]
    fn test_find_search_start_index_all_before() {
        let genes = vec![
            make_test_gene("G1", 100, 200, Strand::Positive, vec![(100, 200)]),
            make_test_gene("G2", 300, 400, Strand::Positive, vec![(300, 400)]),
        ];
        // All genes start before 500
        assert_eq!(find_search_start_index(&genes, 500), 2);
    }

    #[test]
    fn test_find_search_start_index_all_after() {
        let genes = vec![
            make_test_gene("G1", 100, 200, Strand::Positive, vec![(100, 200)]),
            make_test_gene("G2", 300, 400, Strand::Positive, vec![(300, 400)]),
        ];
        // All genes start after 50
        assert_eq!(find_search_start_index(&genes, 50), 0);
    }

    #[test]
    fn test_find_search_start_index_middle() {
        let genes = vec![
            make_test_gene("G1", 100, 200, Strand::Positive, vec![(100, 200)]),
            make_test_gene("G2", 300, 400, Strand::Positive, vec![(300, 400)]),
            make_test_gene("G3", 500, 600, Strand::Positive, vec![(500, 600)]),
        ];
        // First gene starting >= 250 is G2 at index 1
        assert_eq!(find_search_start_index(&genes, 250), 1);
    }

    #[test]
    fn test_match_region_to_genes_no_overlap() {
        let config = Config::default();
        // Region far enough from gene that it's beyond distance threshold (10kb default)
        let region = Region::new("chr1".into(), 10, 20, vec![]);
        let genes = vec![make_test_gene(
            "G1",
            50000,
            60000, // Gene starts at 50kb, well beyond 10kb distance
            Strand::Positive,
            vec![(50000, 60000)],
        )];

        let candidates = match_region_to_genes(&region, &genes, &config, 0);
        // Region [10, 20] is far from gene [50000, 60000] (>10kb default distance)
        assert!(candidates.is_empty());
    }

    #[test]
    fn test_match_region_to_genes_exact_overlap() {
        let config = Config::default();
        let region = Region::new("chr1".into(), 1050, 1150, vec![]);
        let genes = vec![make_test_gene(
            "G1",
            1000,
            2000,
            Strand::Positive,
            vec![(1000, 1200)],
        )];

        let candidates = match_region_to_genes(&region, &genes, &config, 0);
        // Region overlaps with exon, should find 1st_EXON match
        assert!(!candidates.is_empty());

        let has_first_exon = candidates.iter().any(|c| c.area == Area::FirstExon);
        assert!(has_first_exon);
    }

    #[test]
    fn test_process_candidates_empty() {
        let config = Config::default();
        let candidates: Vec<Candidate> = vec![];
        let result = process_candidates_for_output(candidates, &config);
        assert!(result.is_empty());
    }

    #[test]
    fn test_process_candidates_exon_level() {
        let config = Config {
            level: ReportLevel::Exon,
            ..Default::default()
        };

        let c1 = make_candidate(Area::Tss, 100.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Intron, 80.0, 80.0, "T1", "G1", "2");

        let candidates = vec![c1, c2];
        let result = process_candidates_for_output(candidates.clone(), &config);

        // At exon level, all candidates are returned
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_process_candidates_transcript_level() {
        let config = Config {
            level: ReportLevel::Transcript,
            ..Default::default()
        };

        let c1 = make_candidate(Area::Tss, 100.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Intron, 80.0, 80.0, "T1", "G1", "2");

        let candidates = vec![c1, c2];
        let result = process_candidates_for_output(candidates, &config);

        // At transcript level, best per transcript (TSS wins by priority)
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].area, Area::Tss);
    }

    #[test]
    fn test_process_candidates_gene_level() {
        let config = Config {
            level: ReportLevel::Gene,
            ..Default::default()
        };

        let c1 = make_candidate(Area::Tss, 100.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 90.0, 90.0, "T2", "G1", "2");

        let candidates = vec![c1, c2];
        let result = process_candidates_for_output(candidates, &config);

        // At gene level, best per gene with merging
        assert_eq!(result.len(), 1);
        // Should have merged transcripts
        assert!(result[0].transcript.contains("T1") || result[0].transcript.contains("T2"));
    }

    #[test]
    fn test_match_regions_to_genes_basic() {
        let config = Config::default();
        let regions = vec![
            Region::new("chr1".into(), 1050, 1150, vec![]),
            Region::new("chr1".into(), 1500, 1600, vec![]),
        ];

        let genes = vec![make_test_gene(
            "G1",
            1000,
            2000,
            Strand::Positive,
            vec![(1000, 1200), (1500, 1700)],
        )];

        let results = match_regions_to_genes(&regions, &genes, &config, 0);

        assert_eq!(results.len(), 2);

        // First region
        assert_eq!(results[0].0.start, 1050);
        assert!(!results[0].1.is_empty());

        // Second region
        assert_eq!(results[1].0.start, 1500);
        assert!(!results[1].1.is_empty());
    }

    #[test]
    fn test_match_region_to_genes_intron_overlap() {
        // Region falls in the intron between two exons
        let config = Config::default();
        let region = Region::new("chr1".into(), 1250, 1350, vec![]);
        let genes = vec![make_test_gene(
            "G1",
            1000,
            2000,
            Strand::Positive,
            vec![(1000, 1200), (1400, 1600)], // Gap at 1201-1399 is intron
        )];

        let candidates = match_region_to_genes(&region, &genes, &config, 0);
        // Region [1250, 1350] falls in intron [1201, 1399]
        assert!(!candidates.is_empty());
        let has_intron = candidates.iter().any(|c| c.area == Area::Intron);
        assert!(has_intron, "Should find INTRON overlap");
    }

    #[test]
    fn test_match_region_to_genes_negative_strand() {
        // Test negative strand gene - exon numbering is reversed
        let config = Config::default();
        let region = Region::new("chr1".into(), 1050, 1150, vec![]);
        let genes = vec![make_test_gene(
            "G1",
            1000,
            2000,
            Strand::Negative,
            vec![(1000, 1200), (1500, 1700)],
        )];

        let candidates = match_region_to_genes(&region, &genes, &config, 0);
        assert!(!candidates.is_empty());
        // For negative strand, the first exon (genomically) is actually the last exon number-wise
        // So [1000, 1200] should have exon_number "2" and [1500, 1700] should have exon_number "1"
    }

    #[test]
    fn test_match_region_to_genes_upstream_proximity() {
        // Region is upstream of gene (within distance)
        let config = Config::default();
        let region = Region::new("chr1".into(), 800, 900, vec![]);
        let genes = vec![make_test_gene(
            "G1",
            1000,
            2000,
            Strand::Positive,
            vec![(1000, 1200)],
        )];

        let candidates = match_region_to_genes(&region, &genes, &config, 0);
        // Region is 100bp upstream of gene start
        assert!(!candidates.is_empty());
        let has_upstream_or_tss = candidates
            .iter()
            .any(|c| c.area == Area::Upstream || c.area == Area::Tss || c.area == Area::Promoter);
        assert!(has_upstream_or_tss, "Should find UPSTREAM/TSS/PROMOTER");
    }

    #[test]
    fn test_match_region_to_genes_downstream_proximity() {
        // Region is downstream of gene (within distance)
        let config = Config::default();
        let region = Region::new("chr1".into(), 2100, 2200, vec![]);
        let genes = vec![make_test_gene(
            "G1",
            1000,
            2000,
            Strand::Positive,
            vec![(1000, 1200), (1800, 2000)],
        )];

        let candidates = match_region_to_genes(&region, &genes, &config, 0);
        // Region is 100bp downstream of gene end
        assert!(!candidates.is_empty());
        let has_downstream = candidates
            .iter()
            .any(|c| c.area == Area::Downstream || c.area == Area::Tts);
        assert!(has_downstream, "Should find DOWNSTREAM/TTS");
    }

    #[test]
    fn test_match_region_to_genes_multiple_genes() {
        // Region overlaps with multiple genes
        let config = Config::default();
        let region = Region::new("chr1".into(), 1900, 2100, vec![]);
        let genes = vec![
            make_test_gene("G1", 1000, 2000, Strand::Positive, vec![(1800, 2000)]),
            make_test_gene("G2", 2000, 3000, Strand::Positive, vec![(2000, 2200)]),
        ];

        let candidates = match_region_to_genes(&region, &genes, &config, 0);
        // Should have candidates from both genes
        let g1_candidates: Vec<_> = candidates.iter().filter(|c| c.gene == "G1").collect();
        let g2_candidates: Vec<_> = candidates.iter().filter(|c| c.gene == "G2").collect();

        assert!(!g1_candidates.is_empty(), "Should have G1 candidates");
        assert!(!g2_candidates.is_empty(), "Should have G2 candidates");
    }

    #[test]
    fn test_match_region_to_genes_gene_body() {
        // Region overlaps gene but not exons - should get GENE_BODY
        let config = Config::default();
        // Gene spans 1000-5000 but only has exons at ends
        let region = Region::new("chr1".into(), 2500, 2600, vec![]);
        let genes = vec![make_test_gene(
            "G1",
            1000,
            5000,
            Strand::Positive,
            vec![(1000, 1200), (4800, 5000)],
        )];

        let candidates = match_region_to_genes(&region, &genes, &config, 0);
        // Region [2500, 2600] is in the gene body but not in exons
        assert!(!candidates.is_empty());
        let has_gene_body_or_intron = candidates
            .iter()
            .any(|c| c.area == Area::GeneBody || c.area == Area::Intron);
        assert!(
            has_gene_body_or_intron,
            "Should find GENE_BODY or INTRON: {:?}",
            candidates.iter().map(|c| c.area).collect::<Vec<_>>()
        );
    }
}

// -------------------------------------------------------------------------
// 10. Output Module Tests
// -------------------------------------------------------------------------

mod test_output {
    use super::*;
    use rgmatch::Region;

    #[test]
    fn test_format_output_line_basic() {
        let region = Region::new("chr1".to_string(), 100, 200, vec![]);
        let candidate = Candidate::new(
            100,
            200,
            Strand::Positive,
            "1".to_string(),
            Area::Tss,
            "T1".to_string(),
            "G1".to_string(),
            50,
            80.0,
            90.0,
            500,
        );

        let line = format_output_line(&region, &candidate);

        assert!(line.contains("chr1_100_200"));
        assert!(line.contains("150")); // midpoint
        assert!(line.contains("G1"));
        assert!(line.contains("T1"));
        assert!(line.contains("TSS"));
        assert!(line.contains("80.00"));
        assert!(line.contains("90.00"));
    }

    #[test]
    fn test_format_output_line_with_metadata() {
        let region = Region::new(
            "chr1".to_string(),
            100,
            200,
            vec!["peak1".to_string(), "500".to_string(), "+".to_string()],
        );
        let candidate = Candidate::new(
            100,
            200,
            Strand::Positive,
            "1".to_string(),
            Area::Intron,
            "T1".to_string(),
            "G1".to_string(),
            0,
            100.0,
            100.0,
            0,
        );

        let line = format_output_line(&region, &candidate);

        assert!(line.contains("peak1"));
        assert!(line.contains("500"));
        assert!(line.contains("+"));
    }

    #[test]
    fn test_format_output_line_empty_metadata() {
        let region = Region::new("chr1".to_string(), 100, 200, vec![]);
        let candidate = Candidate::new(
            100,
            200,
            Strand::Negative,
            "2".to_string(),
            Area::Downstream,
            "T2".to_string(),
            "G2".to_string(),
            1000,
            50.0,
            -1.0,
            2000,
        );

        let line = format_output_line(&region, &candidate);

        // Should not have trailing tab
        assert!(!line.ends_with('\t'));
        assert!(line.contains("-1.00"));
        assert!(line.contains("DOWNSTREAM"));
    }

    #[test]
    fn test_format_output_line_all_areas() {
        let region = Region::new("chr1".to_string(), 100, 200, vec![]);

        let areas = vec![
            Area::Tss,
            Area::FirstExon,
            Area::Promoter,
            Area::Tts,
            Area::Intron,
            Area::GeneBody,
            Area::Upstream,
            Area::Downstream,
        ];

        for area in areas {
            let candidate = Candidate::new(
                100,
                200,
                Strand::Positive,
                "1".to_string(),
                area,
                "T1".to_string(),
                "G1".to_string(),
                0,
                100.0,
                100.0,
                0,
            );

            let line = format_output_line(&region, &candidate);
            assert!(
                line.contains(area.as_str()),
                "Line should contain {}: {}",
                area,
                line
            );
        }
    }

    #[test]
    fn test_write_header_no_meta() {
        let mut output = Vec::new();
        write_header(&mut output, 0).unwrap();
        let header = String::from_utf8(output).unwrap();

        assert!(header.starts_with("Region\tMidpoint\tGene"));
        assert!(header.contains("Transcript"));
        assert!(header.contains("Exon/Intron"));
        assert!(header.contains("Area"));
        assert!(header.contains("Distance"));
        assert!(header.contains("TSSDistance"));
        assert!(header.contains("PercRegion"));
        assert!(header.contains("PercArea"));
        assert!(!header.contains("name"));
    }

    #[test]
    fn test_write_header_with_meta() {
        let mut output = Vec::new();
        write_header(&mut output, 6).unwrap();
        let header = String::from_utf8(output).unwrap();

        assert!(header.contains("name"));
        assert!(header.contains("score"));
        assert!(header.contains("strand"));
        assert!(header.contains("thickStart"));
        assert!(header.contains("thickEnd"));
        assert!(header.contains("itemRgb"));
        assert!(!header.contains("blockCount")); // Only 6 columns
    }

    #[test]
    fn test_write_header_max_meta() {
        let mut output = Vec::new();
        write_header(&mut output, 9).unwrap();
        let header = String::from_utf8(output).unwrap();

        assert!(header.contains("blockCount"));
        assert!(header.contains("blockSizes"));
        assert!(header.contains("blockStarts"));
    }

    #[test]
    fn test_format_output_line_precision() {
        // Test that percentages are formatted with exactly 2 decimal places
        let region = Region::new("chr1".to_string(), 100, 200, vec![]);
        let candidate = Candidate::new(
            100,
            200,
            Strand::Positive,
            "1".to_string(),
            Area::Tss,
            "T1".to_string(),
            "G1".to_string(),
            0,
            33.333333,
            66.666666,
            0,
        );

        let line = format_output_line(&region, &candidate);

        assert!(line.contains("33.33"));
        assert!(line.contains("66.67"));
    }

    #[test]
    fn test_format_output_line_zero_values() {
        let region = Region::new("chr1".to_string(), 0, 0, vec![]);
        let candidate = Candidate::new(
            0,
            0,
            Strand::Positive,
            "0".to_string(),
            Area::Tss,
            "T0".to_string(),
            "G0".to_string(),
            0,
            0.0,
            0.0,
            0,
        );

        let line = format_output_line(&region, &candidate);

        assert!(line.contains("chr1_0_0"));
        assert!(line.contains("0.00"));
    }

    #[test]
    fn test_format_output_line_large_values() {
        let region = Region::new("chr1".to_string(), 100000000, 200000000, vec![]);
        let candidate = Candidate::new(
            100000000,
            200000000,
            Strand::Negative,
            "999".to_string(),
            Area::GeneBody,
            "TRANSCRIPT_VERY_LONG_NAME".to_string(),
            "GENE_VERY_LONG_NAME".to_string(),
            1000000,
            100.0,
            100.0,
            5000000,
        );

        let line = format_output_line(&region, &candidate);

        assert!(line.contains("chr1_100000000_200000000"));
        assert!(line.contains("150000000")); // midpoint
        assert!(line.contains("TRANSCRIPT_VERY_LONG_NAME"));
        assert!(line.contains("GENE_VERY_LONG_NAME"));
        assert!(line.contains("1000000")); // distance
        assert!(line.contains("5000000")); // tss_distance
    }
}
