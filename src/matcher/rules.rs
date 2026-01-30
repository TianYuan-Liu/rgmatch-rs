//! Rule application logic for selecting best candidates.
//!
//! This module implements the applyRules and selectTranscript functions
//! for filtering and selecting the best candidates based on priority rules.

use ahash::{AHashMap, AHashSet};

use crate::types::{Area, Candidate};

/// Order keys by their first appearance in the candidates list.
///
/// This preserves "insertion order" (file order) to match Python behavior.
/// Keys not found in candidates are sorted and appended at the end.
fn order_keys_by_occurrence<'a, F>(
    candidates: &'a [Candidate],
    grouped_by: &'a AHashMap<String, Vec<usize>>,
    key_fn: F,
) -> Vec<&'a String>
where
    F: Fn(&Candidate) -> &String,
{
    let mut key_order = Vec::new();
    let mut seen = AHashSet::new();

    // Iterate candidates to find unique keys in order of first appearance
    for c in candidates {
        let key = key_fn(c);
        if grouped_by.contains_key(key) && !seen.contains(key) {
            seen.insert(key);
            key_order.push(key);
        }
    }

    // Add any keys from grouped_by that weren't in candidates (unlikely but safe)
    let mut remaining_keys: Vec<&String> =
        grouped_by.keys().filter(|k| !seen.contains(*k)).collect();
    remaining_keys.sort();
    key_order.extend(remaining_keys);

    key_order
}

/// Apply priority rules to select the best candidate per group.
///
/// Filters candidates by percentage thresholds and applies rule-based
/// priority ordering to resolve ties.
///
/// # Arguments
/// * `candidates` - List of Candidate objects to filter
/// * `grouped_by` - Map from group ID to list of candidate indices
/// * `perc_region` - Percentage of region threshold (default 50)
/// * `perc_area` - Percentage of area threshold (default 90)
/// * `rules` - Priority order of areas
///
/// # Returns
/// Filtered list of Candidate objects to report.
pub fn apply_rules(
    candidates: &[Candidate],
    grouped_by: &AHashMap<String, Vec<usize>>,
    perc_region: f64,
    perc_area: f64,
    rules: &[Area],
) -> Vec<Candidate> {
    let mut to_report = Vec::new();

    let key_order = order_keys_by_occurrence(candidates, grouped_by, |c| &c.transcript);

    for key in key_order {
        let positions = &grouped_by[key];
        if positions.len() == 1 {
            to_report.push(candidates[positions[0]].clone());
            continue;
        }

        // Step 1: Filter by %Region threshold
        let mut tmp_results_region: Vec<&Candidate> = positions
            .iter()
            .filter_map(|&pos| {
                let c = &candidates[pos];
                if c.pctg_region >= perc_region {
                    Some(c)
                } else {
                    None
                }
            })
            .collect();

        if tmp_results_region.len() == 1 {
            to_report.push(tmp_results_region[0].clone());
            continue;
        }

        // If none pass, fallback to all candidates
        if tmp_results_region.is_empty() {
            tmp_results_region = positions.iter().map(|&pos| &candidates[pos]).collect();
        }

        if tmp_results_region.len() > 1 {
            // Step 2: Filter by %Area threshold
            let mut tmp_results: Vec<&Candidate> = tmp_results_region
                .iter()
                .filter(|c| c.pctg_area >= perc_area)
                .copied()
                .collect();

            if tmp_results.len() == 1 {
                to_report.push(tmp_results[0].clone());
                continue;
            }

            // If none pass, fallback to all region-filtered candidates
            if tmp_results.is_empty() {
                tmp_results = tmp_results_region;
            }

            if tmp_results.len() > 1 {
                // Step 3: Find max pctg_region among remaining
                let maximum_pctg = tmp_results
                    .iter()
                    .map(|c| c.pctg_region)
                    .fold(0.0_f64, |a, b| a.max(b));

                let region_candidates: Vec<&Candidate> = tmp_results
                    .iter()
                    .filter(|c| c.pctg_region == maximum_pctg)
                    .copied()
                    .collect();

                if region_candidates.len() == 1 {
                    to_report.push(region_candidates[0].clone());
                } else {
                    // Step 4: Apply rules priority order for final selection
                    // Report all that match the first matching rule (ties allowed)
                    let mut found = false;
                    for &area_rule in rules {
                        for &candidate in &region_candidates {
                            if candidate.area == area_rule {
                                to_report.push(candidate.clone());
                                found = true;
                            }
                        }
                        if found {
                            break;
                        }
                    }
                }
            }
        }
    }

    to_report
}

/// Select best transcript from candidates grouped by gene.
///
/// Applies priority rules and merges tied candidates into a single
/// representative with combined transcript/exon information.
///
/// # Arguments
/// * `candidates` - List of Candidate objects to filter
/// * `grouped_by` - Map from gene ID to list of candidate indices
/// * `rules` - Priority order of areas
///
/// # Returns
/// Filtered list of Candidate objects with merged tie information.
pub fn select_transcript(
    candidates: &[Candidate],
    grouped_by: &AHashMap<String, Vec<usize>>,
    rules: &[Area],
) -> Vec<Candidate> {
    let mut to_report = Vec::new();

    // Iterate keys in order of first appearance in candidates (grouped by Gene ID)
    let key_order = order_keys_by_occurrence(candidates, grouped_by, |c| &c.gene);

    for key in key_order {
        let positions = &grouped_by[key];
        if positions.len() == 1 {
            to_report.push(candidates[positions[0]].clone());
            continue;
        }

        // Group by area
        let mut by_area: AHashMap<Area, Vec<usize>> = AHashMap::new();
        for &pos in positions {
            let candidate = &candidates[pos];
            by_area.entry(candidate.area).or_default().push(pos);
        }

        // Apply rules to find winning area
        let mut area_winner: Option<Area> = None;
        for &area_rule in rules {
            if by_area.contains_key(&area_rule) {
                area_winner = Some(area_rule);
                break;
            }
        }

        // Fallback to first available candidate's Area if no rules match
        // "First" means the first one in the list of positions, which preserves order.
        if area_winner.is_none() {
            if let Some(&first_pos) = positions.first() {
                area_winner = Some(candidates[first_pos].area);
            }
        }

        let area_winner = match area_winner {
            Some(a) => a,
            None => continue,
        };

        let winner_positions = &by_area[&area_winner];

        if winner_positions.len() == 1 {
            to_report.push(candidates[winner_positions[0]].clone());
        } else {
            // Merge all tied candidates
            let mut transcripts = String::new();
            let mut exons = String::new();
            let mut max_parea = 0.0_f64;
            let mut max_pregion = 0.0_f64;

            for &pos in winner_positions {
                let c = &candidates[pos];
                transcripts.push_str(&c.transcript);
                transcripts.push(',');
                exons.push_str(&c.exon_number);
                exons.push(',');
                max_parea = max_parea.max(c.pctg_area);
                max_pregion = max_pregion.max(c.pctg_region);
            }

            // Remove trailing comma
            transcripts.pop();
            exons.pop();

            // Use first candidate as reference for other fields
            let ref_candidate = &candidates[winner_positions[0]];
            let merged = Candidate::new(
                ref_candidate.start,
                ref_candidate.end,
                ref_candidate.strand,
                exons,
                ref_candidate.area,
                transcripts,
                ref_candidate.gene.clone(),
                ref_candidate.distance,
                max_pregion,
                max_parea,
                ref_candidate.tss_distance,
            );
            to_report.push(merged);
        }
    }

    to_report
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Strand;

    fn make_candidate(area: Area, pctg_region: f64, pctg_area: f64, transcript: &str) -> Candidate {
        Candidate::new(
            100,
            200,
            Strand::Positive,
            "1".to_string(),
            area,
            transcript.to_string(),
            "G1".to_string(),
            0,
            pctg_region,
            pctg_area,
            100,
        )
    }

    #[test]
    fn test_priority_logic() {
        let rules = vec![
            Area::Tss,
            Area::FirstExon,
            Area::Promoter,
            Area::Tts,
            Area::Intron,
            Area::GeneBody,
            Area::Upstream,
            Area::Downstream,
        ];

        let c1 = make_candidate(Area::Intron, 100.0, 100.0, "T1");
        let c2 = make_candidate(Area::Tss, 100.0, 100.0, "T1"); // Should win
        let c3 = make_candidate(Area::GeneBody, 100.0, 100.0, "T1");

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

        let c1 = make_candidate(Area::Intron, 100.0, 100.0, "T1");
        let c2 = make_candidate(Area::Tss, 100.0, 100.0, "T1");

        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("trans1".to_string(), vec![0, 1]);

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);

        assert_eq!(result.len(), 1);
        assert_eq!(result[0].area, Area::Intron);
    }

    #[test]
    fn test_single_candidate() {
        let rules = vec![Area::Tss];
        let c1 = make_candidate(Area::Tss, 100.0, 100.0, "T1");

        let candidates = vec![c1];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("T1".to_string(), vec![0]);

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);

        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_pctg_region_threshold() {
        let rules = vec![Area::Tss, Area::Intron];

        let c1 = make_candidate(Area::Intron, 60.0, 100.0, "T1"); // Passes
        let c2 = make_candidate(Area::Tss, 40.0, 100.0, "T1"); // Fails threshold

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

        let c1 = make_candidate(Area::Intron, 30.0, 100.0, "T1");
        let c2 = make_candidate(Area::Tss, 40.0, 100.0, "T1");

        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("T1".to_string(), vec![0, 1]);

        let result = apply_rules(&candidates, &grouped_by, 90.0, 90.0, &rules);

        // Should still pick one based on rules priority
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].area, Area::Tss);
    }

    #[test]
    fn test_select_transcript_single() {
        let rules = vec![Area::Tss];
        let c1 = make_candidate(Area::Tss, 100.0, 100.0, "T1");

        let candidates = vec![c1];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("G1".to_string(), vec![0]);

        let result = select_transcript(&candidates, &grouped_by, &rules);

        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_select_transcript_different_areas() {
        let rules = vec![Area::Tss, Area::Intron];

        let c1 = make_candidate(Area::Intron, 100.0, 100.0, "T1");
        let c2 = make_candidate(Area::Tss, 100.0, 100.0, "T2");

        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("G1".to_string(), vec![0, 1]);

        let result = select_transcript(&candidates, &grouped_by, &rules);

        assert_eq!(result.len(), 1);
        assert_eq!(result[0].area, Area::Tss);
    }

    #[test]
    fn test_select_transcript_same_area_merge() {
        let rules = vec![Area::Tss];

        let mut c1 = make_candidate(Area::Tss, 80.0, 70.0, "T1");
        c1.exon_number = "1".to_string();
        let mut c2 = make_candidate(Area::Tss, 90.0, 60.0, "T2");
        c2.exon_number = "2".to_string();

        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("G1".to_string(), vec![0, 1]);

        let result = select_transcript(&candidates, &grouped_by, &rules);

        assert_eq!(result.len(), 1);
        assert!(result[0].transcript.contains("T1"));
        assert!(result[0].transcript.contains("T2"));
        assert!(result[0].exon_number.contains("1"));
        assert!(result[0].exon_number.contains("2"));
        assert_eq!(result[0].pctg_region, 90.0); // max of 80, 90
        assert_eq!(result[0].pctg_area, 70.0); // max of 70, 60
    }

    #[test]
    fn test_max_pctg_region_tiebreaker() {
        let rules = vec![Area::Tss];

        let c1 = make_candidate(Area::Tss, 80.0, 100.0, "T1");
        let c2 = make_candidate(Area::Tss, 90.0, 100.0, "T2"); // Higher

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

        let c1 = make_candidate(Area::Tss, 80.0, 100.0, "T1");
        let c2 = make_candidate(Area::Tss, 80.0, 100.0, "T2");

        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("T1".to_string(), vec![0, 1]);

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);

        // Both should be reported (tie)
        assert_eq!(result.len(), 2);
    }
}
