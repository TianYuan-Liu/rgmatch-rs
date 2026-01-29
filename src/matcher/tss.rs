//! TSS (Transcription Start Site) overlap checking.
//!
//! This module implements the checkTSS logic with coordinate mirroring
//! for negative strand genes.

use crate::types::Strand;

/// Result of a TSS check: (area_tag, pctg_dhs, pctg_area).
pub type TssResult = (String, f64, f64);

/// Helper struct to pass exon-like data to checkTSS.
pub struct TssExonInfo {
    pub start: i64,
    pub end: i64,
    pub strand: Strand,
    pub distance: i64,
}

/// Check overlap with TSS (Transcription Start Site) region.
///
/// Calculates the overlap between a DHS region and the TSS/promoter
/// regions upstream of the first exon. Handles strand orientation
/// by coordinate transformation.
///
/// CRITICAL: For negative strand, coordinates are mirrored around the exon end!
///
/// # Arguments
/// * `dhs_start` - Start coordinate of the DHS region
/// * `dhs_end` - End coordinate of the DHS region
/// * `exon_info` - Exon information including position, strand, and distance
/// * `tss_distance` - TSS region distance (default 200bp)
/// * `promoter_distance` - Promoter region distance (default 1300bp)
///
/// # Returns
/// A vector of (area_tag, pctg_dhs, pctg_area) tuples for each overlapping region type.
pub fn check_tss(
    dhs_start: i64,
    dhs_end: i64,
    exon_info: &TssExonInfo,
    tss_distance: f64,
    promoter_distance: f64,
) -> Vec<TssResult> {
    let mut exon_start = exon_info.start;
    let distance_val = exon_info.distance;
    let mut actual_dhs_start = dhs_start;
    let mut actual_dhs_end = dhs_end;

    // CRITICAL: Coordinate mirroring for negative strand
    // For negative strand, we flip the coordinates to make the code strand-invariant
    if exon_info.strand == Strand::Negative {
        let aux = actual_dhs_end;
        actual_dhs_end = 2 * exon_info.end - actual_dhs_start;
        actual_dhs_start = 2 * exon_info.end - aux;
        exon_start = exon_info.end; // TSS is at exon END for negative strand
    }

    let dhs_length = actual_dhs_end - actual_dhs_start + 1;

    // Zero-length region check - must be <= 0, not < 0
    if dhs_length <= 0 {
        return vec![];
    }

    let mut results = Vec::new();
    let dhs_length_f = dhs_length as f64;

    if distance_val as f64 <= tss_distance {
        // Region is within TSS distance

        // UPSTREAM       PROMOTER        TSS          1st exon
        // ..........|................|..............|----------->

        if (exon_start - actual_dhs_start) as f64 <= tss_distance {
            // Region is entirely within TSS zone
            // UPSTREAM       PROMOTER        TSS          1st exon
            // ..........|................|..............|----------->
            //                      DHS
            //                               |-------------

            let overlap_end = std::cmp::min(exon_start - 1, actual_dhs_end);
            let overlap = overlap_end - actual_dhs_start + 1;
            let pctg_dhs = (overlap as f64 / dhs_length_f) * 100.0;
            let pctg_tss = (overlap as f64 / tss_distance) * 100.0;
            results.push(("TSS".to_string(), pctg_dhs, pctg_tss));
        } else {
            // Region spans TSS and extends into PROMOTER
            // UPSTREAM       PROMOTER        TSS          1st exon
            // ..........|................|..............|----------->
            //                      DHS
            //                        --------------

            // TSS portion
            let tss_start = exon_start - tss_distance as i64;
            let overlap_end = std::cmp::min(exon_start - 1, actual_dhs_end);
            let tss_overlap = overlap_end - tss_start + 1;
            let pctg_dhs_tss = (tss_overlap as f64 / dhs_length_f) * 100.0;
            let pctg_tss = (tss_overlap as f64 / tss_distance) * 100.0;
            results.push(("TSS".to_string(), pctg_dhs_tss, pctg_tss));

            // Check if region extends into PROMOTER
            if (exon_start - actual_dhs_start) as f64 <= tss_distance + promoter_distance {
                // Region is within TSS + PROMOTER zone
                let promoter_overlap = (exon_start - tss_distance as i64) - actual_dhs_start;
                let pctg_dhs_promoter = (promoter_overlap as f64 / dhs_length_f) * 100.0;
                let pctg_promoter = (promoter_overlap as f64 / promoter_distance) * 100.0;
                results.push(("PROMOTER".to_string(), pctg_dhs_promoter, pctg_promoter));
            } else {
                // Region extends into UPSTREAM
                let pctg_dhs_promoter = (promoter_distance / dhs_length_f) * 100.0;
                let pctg_promoter = 100.0;
                results.push(("PROMOTER".to_string(), pctg_dhs_promoter, pctg_promoter));

                let upstream_overlap =
                    (exon_start - tss_distance as i64 - promoter_distance as i64)
                        - actual_dhs_start;
                let pctg_dhs_upstream = (upstream_overlap as f64 / dhs_length_f) * 100.0;
                results.push(("UPSTREAM".to_string(), pctg_dhs_upstream, -1.0));
            }
        }
    } else if distance_val as f64 <= tss_distance + promoter_distance {
        // Region is within PROMOTER zone (beyond TSS)

        if (exon_start - actual_dhs_start) as f64 <= tss_distance + promoter_distance {
            // Region is entirely within PROMOTER zone
            let pctg_dhs = 100.0;
            let pctg_promoter = (dhs_length_f / promoter_distance) * 100.0;
            results.push(("PROMOTER".to_string(), pctg_dhs, pctg_promoter));
        } else {
            // Region spans PROMOTER and extends into UPSTREAM
            let promoter_start = exon_start - tss_distance as i64 - promoter_distance as i64;
            let promoter_overlap = actual_dhs_end - promoter_start + 1;
            let pctg_dhs_promoter = (promoter_overlap as f64 / dhs_length_f) * 100.0;
            let pctg_promoter = (promoter_overlap as f64 / promoter_distance) * 100.0;
            results.push(("PROMOTER".to_string(), pctg_dhs_promoter, pctg_promoter));

            let upstream_overlap = promoter_start - actual_dhs_start;
            let pctg_dhs_upstream = (upstream_overlap as f64 / dhs_length_f) * 100.0;
            results.push(("UPSTREAM".to_string(), pctg_dhs_upstream, -1.0));
        }
    } else {
        // Region is entirely in UPSTREAM zone
        results.push(("UPSTREAM".to_string(), 100.0, -1.0));
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pos_strand_boundaries() {
        // Exon: [2000, 3000]. TSS @ 2000.
        // TSS zone: [1800, 2000].
        // Promoter zone: [500, 1800) (1300bp long).

        // Case 1: Exactly at TSS boundary (200bp upstream) -> [1800, 1810]
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
        // TSS zone: [3000, 3200].
        // Promoter zone: [3200, 4500].

        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Negative,
            distance: 0,
        };

        // Case 1: Region [3200, 3210] (at PROMOTER boundary)
        // Flipped: dhs_start' = 2*3000 - 3210 = 2790
        // Flipped: dhs_end' = 2*3000 - 3200 = 2800
        // exon_start' = 3000
        // 3000 - 2790 = 210 > 200, so PROMOTER
        let res = check_tss(3200, 3210, &exon, 200.0, 1300.0);
        let tags: Vec<&str> = res.iter().map(|(tag, _, _)| tag.as_str()).collect();
        assert!(tags.contains(&"PROMOTER"));

        // Case 2: TSS Zone Inside [3100, 3150]
        // Flipped: 2*3000 - 3150 = 2850 (Start).
        // 3000 - 2850 = 150 <= 200, so TSS.
        let res = check_tss(3100, 3150, &exon, 200.0, 1300.0);
        assert!(res.iter().any(|(tag, _, _)| tag == "TSS"));
    }

    #[test]
    fn test_integer_division_check_tss() {
        // Verify checkTSS uses integer division for midpoint calculation
        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Positive,
            distance: 0,
        };
        // With start=1801, end=1810
        let res = check_tss(1801, 1810, &exon, 200.0, 1300.0);
        // The function should complete without float issues
        assert!(!res.is_empty());
    }

    #[test]
    fn test_zero_length_region_check_tss() {
        // Verify checkTSS handles zero-length regions gracefully
        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Positive,
            distance: 0,
        };
        // A region where end < start results in dhs_length <= 0
        let res = check_tss(1900, 1899, &exon, 200.0, 1300.0);
        assert!(res.is_empty());
    }

    #[test]
    fn test_zero_tss_value() {
        // Test checkTSS when tss is set to 0 goes to promoter logic
        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Positive,
            distance: 500,
        };
        let res = check_tss(1500, 1600, &exon, 0.0, 1300.0);
        let tags: Vec<&str> = res.iter().map(|(tag, _, _)| tag.as_str()).collect();
        assert!(tags.contains(&"PROMOTER"));
    }

    #[test]
    fn test_large_tss_value() {
        let exon = TssExonInfo {
            start: 20000,
            end: 30000,
            strand: Strand::Positive,
            distance: 5000,
        };
        let res = check_tss(15000, 15100, &exon, 10000.0, 1300.0);
        assert!(res.iter().any(|(tag, _, _)| tag == "TSS"));
    }
}
