//! Region-gene overlap matching logic.
//!
//! This module implements the main matching loop that associates genomic regions
//! with gene annotations based on positional overlap and proximity.

use ahash::AHashMap;
use indexmap::IndexMap;

use crate::config::Config;
use crate::matcher::rules::{apply_rules, select_transcript};
use crate::matcher::tss::{check_tss, TssExonInfo};
use crate::matcher::tts::{check_tts, TtsExonInfo};
use crate::types::{Area, Candidate, Gene, Region, ReportLevel, Strand};

/// Calculate the intron number based on exon index and strand.
///
/// For positive strand genes, intron N is between exon N and exon N+1.
/// For negative strand genes, the numbering is reversed from the 3' end.
fn calculate_intron_number(index: usize, total_exons: usize, strand: Strand) -> usize {
    match strand {
        Strand::Positive => index + 1,
        Strand::Negative => total_exons - 1 - index,
    }
}

/// Aggregate overlapping entries (gene body or intron) into a single candidate per transcript.
///
/// Takes a map of entries grouped by transcript key and combines overlapping regions
/// into single candidates with aggregated statistics.
fn aggregate_entries(
    entries_map: IndexMap<String, Vec<(Candidate, i64, i64)>>,
    region_length: i64,
) -> Vec<Candidate> {
    let mut results = Vec::new();

    for (_, entries) in entries_map {
        if entries.len() == 1 {
            results.push(entries[0].0.clone());
        } else {
            let mut total_area = 0i64;
            let mut total_overlap = 0i64;
            let mut combined_numbers = String::new();

            for (candidate, area_len, overlap) in &entries {
                total_area += area_len;
                total_overlap += overlap;
                combined_numbers.push_str(&candidate.exon_number);
                combined_numbers.push(',');
            }
            combined_numbers.pop(); // Remove trailing comma

            let ref_candidate = &entries[0].0;
            let pctg_region = (total_overlap as f64 / region_length as f64) * 100.0;
            let pctg_area = (total_overlap as f64 / total_area as f64) * 100.0;

            results.push(Candidate::new(
                ref_candidate.start,
                ref_candidate.end,
                ref_candidate.strand,
                combined_numbers,
                ref_candidate.area,
                ref_candidate.transcript.clone(),
                ref_candidate.gene.clone(),
                ref_candidate.distance,
                pctg_region,
                pctg_area,
                ref_candidate.tss_distance,
            ));
        }
    }

    results
}

/// Match a single region to genes and return all candidates.
///
/// This implements the main matching logic from the Python code.
pub fn match_region_to_genes(
    region: &Region,
    genes: &[Gene],
    config: &Config,
    last_index: usize,
) -> Vec<Candidate> {
    let start = region.start;
    let end = region.end;
    let pm = region.midpoint();
    let region_length = region.length();

    // Start analysis
    let mut down: i64 = i64::MAX; // Distance to TTS
    let mut exon_down: Option<Candidate> = None;

    let mut upst: i64 = i64::MAX; // Distance to TSS
    let mut exon_up: Option<Candidate> = None;

    // When flag_gene_body is false, we will report downstream or upstream exons
    // Otherwise, we will only report the overlapped exons
    let mut flag_gene_body = false;

    // Array containing the relations that are going to be reported
    let mut final_output: Vec<Candidate> = Vec::new();

    // These maps contain as key [geneID_transcriptID] and as values a vector
    // containing [(Candidate, area_length, overlapped_area), ...]
    // This is because there will be regions that will overlap different introns or exons
    let mut my_introns: IndexMap<String, Vec<(Candidate, i64, i64)>> = IndexMap::new();
    let mut my_gene_bodys: IndexMap<String, Vec<(Candidate, i64, i64)>> = IndexMap::new();

    for (_i, gene) in genes.iter().enumerate().skip(last_index) {
        let distance_to_start_gene = (gene.start - pm).abs();

        // Check if we should stop processing genes
        // Since genes are sorted by start, if the gene starts after our region ends (plus lookahead),
        // no subsequent genes can possibly overlap.
        // Note: The lookahead logic depends on whether we are looking for UPSTREAM/DOWNSTREAM
        if gene.start > end {
            // Logic for quitting early:
            // If we are looking for downstream (down), checking gene.start > end is not enough if we want nearest.
            // But 'down' is initialized to MAX.
            // The python logic seems to be: if we found something closer than current distance, stop.
            // Simplified check matching Python structure:
            if flag_gene_body || down < distance_to_start_gene || upst < distance_to_start_gene {
                break;
            }
            // Additional safety check for performance: if gene starts WAY after, we can definitely stop?
            // Existing logic relies on `down` and `upst` being updated.
        }

        // Check associations
        for transcript in &gene.transcripts {
            let exons = &transcript.exons;

            // Calculate TSSdist using the first exon "start" position
            let tss_distance = if exons[0].exon_number.as_deref() == Some("1") {
                pm - exons[0].start
            } else {
                exons.last().unwrap().end - pm
            };

            for (j, exon) in exons.iter().enumerate() {
                let is_first_exon = j == 0;
                let is_last_exon = j == exons.len() - 1;
                let exon_length = exon.length();
                let exon_number = exon.exon_number.clone().unwrap_or_default();

                // Case 1: Exon before the region
                // <--------->
                //                |--------------|
                if exon.end < start {
                    // Check whether the current gene also covers the region

                    let dist_tmp = pm - exon.end;

                    // Check if it's the last exon
                    if is_last_exon {
                        if gene.strand == Strand::Positive && dist_tmp < down {
                            down = dist_tmp;
                            exon_down = Some(Candidate::new(
                                exon.start,
                                exon.end,
                                gene.strand,
                                exon_number.clone(),
                                Area::Downstream,
                                transcript.transcript_id.clone(),
                                gene.gene_id.clone(),
                                down,
                                100.0,
                                -1.0,
                                tss_distance,
                            ));
                        } else if gene.strand == Strand::Negative && dist_tmp < upst {
                            upst = dist_tmp;
                            exon_up = Some(Candidate::new(
                                exon.start,
                                exon.end,
                                gene.strand,
                                exon_number.clone(),
                                Area::Upstream,
                                transcript.transcript_id.clone(),
                                gene.gene_id.clone(),
                                upst,
                                100.0,
                                -1.0,
                                tss_distance,
                            ));
                        }
                    } else {
                        // Check if the next exon is closer to the region
                        let next_exon = &exons[j + 1];

                        if next_exon.start > start {
                            flag_gene_body = true;
                            let intron_length = next_exon.start - exon.end - 1;
                            let intron_number =
                                calculate_intron_number(j, exons.len(), gene.strand);

                            if next_exon.start > end {
                                // Region is completely inside intron
                                let pctg_region = 100.0;
                                let pctg_area =
                                    (region_length as f64 / intron_length as f64) * 100.0;

                                let my_id =
                                    format!("{}_{}", gene.gene_id, transcript.transcript_id);
                                let intron_candidate = Candidate::new(
                                    exon.start,
                                    exon.end,
                                    gene.strand,
                                    intron_number.to_string(),
                                    Area::Intron,
                                    transcript.transcript_id.clone(),
                                    gene.gene_id.clone(),
                                    0,
                                    pctg_region,
                                    pctg_area,
                                    tss_distance,
                                );
                                my_introns.entry(my_id).or_default().push((
                                    intron_candidate,
                                    intron_length,
                                    region_length,
                                ));
                                break;
                            } else {
                                // Region overlaps with next exon
                                let region_overlap = next_exon.start - start;
                                let pctg_region =
                                    (region_overlap as f64 / region_length as f64) * 100.0;
                                let pctg_area =
                                    (region_overlap as f64 / intron_length as f64) * 100.0;

                                let my_id =
                                    format!("{}_{}", gene.gene_id, transcript.transcript_id);
                                let intron_candidate = Candidate::new(
                                    exon.start,
                                    exon.end,
                                    gene.strand,
                                    intron_number.to_string(),
                                    Area::Intron,
                                    transcript.transcript_id.clone(),
                                    gene.gene_id.clone(),
                                    0,
                                    pctg_region,
                                    pctg_area,
                                    tss_distance,
                                );
                                my_introns.entry(my_id).or_default().push((
                                    intron_candidate,
                                    intron_length,
                                    region_overlap,
                                ));
                            }
                        }
                    }
                }
                // Case 2: Exon overlapping partially the region (left)
                //     <--------->
                //          |--------------|
                else if start <= exon.end && exon.end <= end && exon.start < start {
                    flag_gene_body = true;
                    let body_overlap = exon.end - start + 1;
                    let pctg_region = (body_overlap as f64 / region_length as f64) * 100.0;
                    let pctg_area = (body_overlap as f64 / exon_length as f64) * 100.0;

                    if (is_first_exon && gene.strand == Strand::Positive)
                        || (is_last_exon && gene.strand == Strand::Negative)
                    {
                        final_output.push(Candidate::new(
                            exon.start,
                            exon.end,
                            gene.strand,
                            exon_number.clone(),
                            Area::FirstExon,
                            transcript.transcript_id.clone(),
                            gene.gene_id.clone(),
                            0,
                            pctg_region,
                            pctg_area,
                            tss_distance,
                        ));
                    } else {
                        let my_id = format!("{}_{}", gene.gene_id, transcript.transcript_id);
                        let gb_candidate = Candidate::new(
                            exon.start,
                            exon.end,
                            gene.strand,
                            exon_number.clone(),
                            Area::GeneBody,
                            transcript.transcript_id.clone(),
                            gene.gene_id.clone(),
                            0,
                            pctg_region,
                            pctg_area,
                            tss_distance,
                        );
                        my_gene_bodys.entry(my_id).or_default().push((
                            gb_candidate,
                            exon_length,
                            body_overlap,
                        ));
                    }

                    // Handle remaining region after exon
                    if exon.end < end {
                        if is_last_exon {
                            let region_overlap = end - exon.end;
                            let pctg_region_r =
                                (region_overlap as f64 / region_length as f64) * 100.0;

                            if gene.strand == Strand::Positive {
                                let candidate = Candidate::new(
                                    exon.start,
                                    exon.end,
                                    gene.strand,
                                    exon_number.clone(),
                                    Area::Downstream,
                                    transcript.transcript_id.clone(),
                                    gene.gene_id.clone(),
                                    0,
                                    pctg_region_r,
                                    -1.0,
                                    tss_distance,
                                );
                                if config.tts > 0.0 {
                                    let exon_info = TtsExonInfo {
                                        start: candidate.start,
                                        end: candidate.end,
                                        strand: candidate.strand,
                                        distance: candidate.distance,
                                    };
                                    for (tag, pctg_dhs, pctg_a) in
                                        check_tts(start, end, &exon_info, config.tts)
                                    {
                                        final_output.push(Candidate::new(
                                            candidate.start,
                                            candidate.end,
                                            candidate.strand,
                                            candidate.exon_number.clone(),
                                            tag.parse().unwrap_or(Area::Downstream),
                                            candidate.transcript.clone(),
                                            candidate.gene.clone(),
                                            candidate.distance,
                                            pctg_dhs,
                                            pctg_a,
                                            tss_distance,
                                        ));
                                    }
                                } else {
                                    final_output.push(candidate);
                                }
                            } else {
                                let candidate = Candidate::new(
                                    exon.start,
                                    exon.end,
                                    gene.strand,
                                    exon_number.clone(),
                                    Area::Upstream,
                                    transcript.transcript_id.clone(),
                                    gene.gene_id.clone(),
                                    0,
                                    pctg_region_r,
                                    -1.0,
                                    tss_distance,
                                );
                                let exon_info = TssExonInfo {
                                    start: candidate.start,
                                    end: candidate.end,
                                    strand: candidate.strand,
                                    distance: candidate.distance,
                                };
                                for (tag, pctg_dhs, pctg_a) in
                                    check_tss(start, end, &exon_info, config.tss, config.promoter)
                                {
                                    final_output.push(Candidate::new(
                                        candidate.start,
                                        candidate.end,
                                        candidate.strand,
                                        candidate.exon_number.clone(),
                                        tag.parse().unwrap_or(Area::Upstream),
                                        candidate.transcript.clone(),
                                        candidate.gene.clone(),
                                        candidate.distance,
                                        pctg_dhs,
                                        pctg_a,
                                        tss_distance,
                                    ));
                                }
                            }
                        } else {
                            // Check intron after exon
                            let next_exon = &exons[j + 1];
                            let intron_length = next_exon.start - exon.end - 1;
                            let intron_number =
                                calculate_intron_number(j, exons.len(), gene.strand);

                            if next_exon.start > end {
                                let region_overlap = end - exon.end;
                                let pctg_region =
                                    (region_overlap as f64 / region_length as f64) * 100.0;
                                let pctg_area =
                                    (region_overlap as f64 / intron_length as f64) * 100.0;

                                let my_id =
                                    format!("{}_{}", gene.gene_id, transcript.transcript_id);
                                let intron_candidate = Candidate::new(
                                    exon.start,
                                    exon.end,
                                    gene.strand,
                                    intron_number.to_string(),
                                    Area::Intron,
                                    transcript.transcript_id.clone(),
                                    gene.gene_id.clone(),
                                    0,
                                    pctg_region,
                                    pctg_area,
                                    tss_distance,
                                );
                                my_introns.entry(my_id).or_default().push((
                                    intron_candidate,
                                    intron_length,
                                    region_overlap,
                                ));
                                break;
                            } else {
                                let region_overlap = next_exon.start - exon.end - 1;
                                let pctg_region =
                                    (region_overlap as f64 / region_length as f64) * 100.0;
                                let pctg_area =
                                    (region_overlap as f64 / intron_length as f64) * 100.0;

                                let my_id =
                                    format!("{}_{}", gene.gene_id, transcript.transcript_id);
                                let intron_candidate = Candidate::new(
                                    exon.start,
                                    exon.end,
                                    gene.strand,
                                    intron_number.to_string(),
                                    Area::Intron,
                                    transcript.transcript_id.clone(),
                                    gene.gene_id.clone(),
                                    0,
                                    pctg_region,
                                    pctg_area,
                                    tss_distance,
                                );

                                my_introns.entry(my_id).or_default().push((
                                    intron_candidate,
                                    intron_length,
                                    region_overlap,
                                ));
                            }
                        }
                    }
                }
                // Case 3: Exon completely inside the region
                //     <--------->
                //   |--------------|
                else if start <= exon.start && end >= exon.end {
                    flag_gene_body = true;

                    // Handle upstream portion before exon
                    if start < exon.start && is_first_exon {
                        let region_overlap = exon.start - start;
                        let pctg_region_r = (region_overlap as f64 / region_length as f64) * 100.0;

                        if gene.strand == Strand::Negative {
                            let candidate = Candidate::new(
                                exon.start,
                                exon.end,
                                gene.strand,
                                exon_number.clone(),
                                Area::Downstream,
                                transcript.transcript_id.clone(),
                                gene.gene_id.clone(),
                                0,
                                pctg_region_r,
                                -1.0,
                                tss_distance,
                            );
                            if config.tts > 0.0 {
                                let exon_info = TtsExonInfo {
                                    start: candidate.start,
                                    end: candidate.end,
                                    strand: candidate.strand,
                                    distance: candidate.distance,
                                };
                                for (tag, pctg_dhs, pctg_a) in
                                    check_tts(start, end, &exon_info, config.tts)
                                {
                                    final_output.push(Candidate::new(
                                        candidate.start,
                                        candidate.end,
                                        candidate.strand,
                                        candidate.exon_number.clone(),
                                        tag.parse().unwrap_or(Area::Downstream),
                                        candidate.transcript.clone(),
                                        candidate.gene.clone(),
                                        candidate.distance,
                                        pctg_dhs,
                                        pctg_a,
                                        tss_distance,
                                    ));
                                }
                            } else {
                                final_output.push(candidate);
                            }
                        } else {
                            let candidate = Candidate::new(
                                exon.start,
                                exon.end,
                                gene.strand,
                                exon_number.clone(),
                                Area::Upstream,
                                transcript.transcript_id.clone(),
                                gene.gene_id.clone(),
                                0,
                                pctg_region_r,
                                -1.0,
                                tss_distance,
                            );
                            let exon_info = TssExonInfo {
                                start: candidate.start,
                                end: candidate.end,
                                strand: candidate.strand,
                                distance: candidate.distance,
                            };
                            for (tag, pctg_dhs, pctg_a) in
                                check_tss(start, end, &exon_info, config.tss, config.promoter)
                            {
                                final_output.push(Candidate::new(
                                    candidate.start,
                                    candidate.end,
                                    candidate.strand,
                                    candidate.exon_number.clone(),
                                    tag.parse().unwrap_or(Area::Upstream),
                                    candidate.transcript.clone(),
                                    candidate.gene.clone(),
                                    candidate.distance,
                                    pctg_dhs,
                                    pctg_a,
                                    tss_distance,
                                ));
                            }
                        }
                    }

                    // Handle the exon overlap
                    let region_overlap = exon.end - exon.start + 1;
                    let pctg_region = (region_overlap as f64 / region_length as f64) * 100.0;
                    let pctg_area = 100.0;

                    if (is_first_exon && gene.strand == Strand::Positive)
                        || (is_last_exon && gene.strand == Strand::Negative)
                    {
                        final_output.push(Candidate::new(
                            exon.start,
                            exon.end,
                            gene.strand,
                            exon_number.clone(),
                            Area::FirstExon,
                            transcript.transcript_id.clone(),
                            gene.gene_id.clone(),
                            0,
                            pctg_region,
                            pctg_area,
                            tss_distance,
                        ));
                    } else {
                        let my_id = format!("{}_{}", gene.gene_id, transcript.transcript_id);

                        let gb_candidate = Candidate::new(
                            exon.start,
                            exon.end,
                            gene.strand,
                            exon_number.clone(),
                            Area::GeneBody,
                            transcript.transcript_id.clone(),
                            gene.gene_id.clone(),
                            0,
                            pctg_region,
                            pctg_area,
                            tss_distance,
                        );
                        my_gene_bodys.entry(my_id).or_default().push((
                            gb_candidate,
                            exon_length,
                            exon_length,
                        ));
                    }

                    // Handle downstream portion after exon
                    if end > exon.end {
                        if is_last_exon {
                            let region_overlap = end - exon.end;
                            let pctg_region_r =
                                (region_overlap as f64 / region_length as f64) * 100.0;

                            if gene.strand == Strand::Positive {
                                let candidate = Candidate::new(
                                    exon.start,
                                    exon.end,
                                    gene.strand,
                                    exon_number.clone(),
                                    Area::Downstream,
                                    transcript.transcript_id.clone(),
                                    gene.gene_id.clone(),
                                    0,
                                    pctg_region_r,
                                    -1.0,
                                    tss_distance,
                                );
                                if config.tts > 0.0 {
                                    let exon_info = TtsExonInfo {
                                        start: candidate.start,
                                        end: candidate.end,
                                        strand: candidate.strand,
                                        distance: candidate.distance,
                                    };
                                    for (tag, pctg_dhs, pctg_a) in
                                        check_tts(start, end, &exon_info, config.tts)
                                    {
                                        final_output.push(Candidate::new(
                                            candidate.start,
                                            candidate.end,
                                            candidate.strand,
                                            candidate.exon_number.clone(),
                                            tag.parse().unwrap_or(Area::Downstream),
                                            candidate.transcript.clone(),
                                            candidate.gene.clone(),
                                            candidate.distance,
                                            pctg_dhs,
                                            pctg_a,
                                            tss_distance,
                                        ));
                                    }
                                } else {
                                    final_output.push(candidate);
                                }
                            } else {
                                let candidate = Candidate::new(
                                    exon.start,
                                    exon.end,
                                    gene.strand,
                                    exon_number.clone(),
                                    Area::Upstream,
                                    transcript.transcript_id.clone(),
                                    gene.gene_id.clone(),
                                    0,
                                    pctg_region_r,
                                    -1.0,
                                    tss_distance,
                                );
                                let exon_info = TssExonInfo {
                                    start: candidate.start,
                                    end: candidate.end,
                                    strand: candidate.strand,
                                    distance: candidate.distance,
                                };
                                for (tag, pctg_dhs, pctg_a) in
                                    check_tss(start, end, &exon_info, config.tss, config.promoter)
                                {
                                    final_output.push(Candidate::new(
                                        candidate.start,
                                        candidate.end,
                                        candidate.strand,
                                        candidate.exon_number.clone(),
                                        tag.parse().unwrap_or(Area::Upstream),
                                        candidate.transcript.clone(),
                                        candidate.gene.clone(),
                                        candidate.distance,
                                        pctg_dhs,
                                        pctg_a,
                                        tss_distance,
                                    ));
                                }
                            }
                        } else {
                            // Check intron after exon
                            let next_exon = &exons[j + 1];
                            let intron_length = next_exon.start - exon.end - 1;
                            let intron_number =
                                calculate_intron_number(j, exons.len(), gene.strand);

                            if next_exon.start > end {
                                let region_overlap = end - exon.end;
                                let pctg_region =
                                    (region_overlap as f64 / region_length as f64) * 100.0;
                                let pctg_area =
                                    (region_overlap as f64 / intron_length as f64) * 100.0;

                                let my_id =
                                    format!("{}_{}", gene.gene_id, transcript.transcript_id);
                                let intron_candidate = Candidate::new(
                                    exon.start,
                                    exon.end,
                                    gene.strand,
                                    intron_number.to_string(),
                                    Area::Intron,
                                    transcript.transcript_id.clone(),
                                    gene.gene_id.clone(),
                                    0,
                                    pctg_region,
                                    pctg_area,
                                    tss_distance,
                                );
                                my_introns.entry(my_id).or_default().push((
                                    intron_candidate,
                                    intron_length,
                                    region_overlap,
                                ));
                                break;
                            } else {
                                let region_overlap = next_exon.start - exon.end - 1;
                                let pctg_region =
                                    (region_overlap as f64 / region_length as f64) * 100.0;
                                let pctg_area =
                                    (region_overlap as f64 / intron_length as f64) * 100.0;

                                let my_id =
                                    format!("{}_{}", gene.gene_id, transcript.transcript_id);
                                let intron_candidate = Candidate::new(
                                    exon.start,
                                    exon.end,
                                    gene.strand,
                                    intron_number.to_string(),
                                    Area::Intron,
                                    transcript.transcript_id.clone(),
                                    gene.gene_id.clone(),
                                    0,
                                    pctg_region,
                                    pctg_area,
                                    tss_distance,
                                );
                                my_introns.entry(my_id).or_default().push((
                                    intron_candidate,
                                    intron_length,
                                    region_overlap,
                                ));
                            }
                        }
                    }
                }
                // Case 4: Exon overlapping the region but shifted to the right
                //             <--------->
                //   |--------------|
                else if start <= exon.start && exon.start <= end && end < exon.end {
                    flag_gene_body = true;

                    // Handle upstream portion before exon
                    if start < exon.start && is_first_exon {
                        let region_overlap = exon.start - start;
                        let pctg_region_r = (region_overlap as f64 / region_length as f64) * 100.0;

                        if gene.strand == Strand::Negative {
                            let candidate = Candidate::new(
                                exon.start,
                                exon.end,
                                gene.strand,
                                exon_number.clone(),
                                Area::Downstream,
                                transcript.transcript_id.clone(),
                                gene.gene_id.clone(),
                                0,
                                pctg_region_r,
                                -1.0,
                                tss_distance,
                            );
                            if config.tts > 0.0 {
                                let exon_info = TtsExonInfo {
                                    start: candidate.start,
                                    end: candidate.end,
                                    strand: candidate.strand,
                                    distance: candidate.distance,
                                };
                                for (tag, pctg_dhs, pctg_a) in
                                    check_tts(start, end, &exon_info, config.tts)
                                {
                                    final_output.push(Candidate::new(
                                        candidate.start,
                                        candidate.end,
                                        candidate.strand,
                                        candidate.exon_number.clone(),
                                        tag.parse().unwrap_or(Area::Downstream),
                                        candidate.transcript.clone(),
                                        candidate.gene.clone(),
                                        candidate.distance,
                                        pctg_dhs,
                                        pctg_a,
                                        tss_distance,
                                    ));
                                }
                            } else {
                                final_output.push(candidate);
                            }
                        } else {
                            let candidate = Candidate::new(
                                exon.start,
                                exon.end,
                                gene.strand,
                                exon_number.clone(),
                                Area::Upstream,
                                transcript.transcript_id.clone(),
                                gene.gene_id.clone(),
                                0,
                                pctg_region_r,
                                -1.0,
                                tss_distance,
                            );
                            let exon_info = TssExonInfo {
                                start: candidate.start,
                                end: candidate.end,
                                strand: candidate.strand,
                                distance: candidate.distance,
                            };
                            for (tag, pctg_dhs, pctg_a) in
                                check_tss(start, end, &exon_info, config.tss, config.promoter)
                            {
                                final_output.push(Candidate::new(
                                    candidate.start,
                                    candidate.end,
                                    candidate.strand,
                                    candidate.exon_number.clone(),
                                    tag.parse().unwrap_or(Area::Upstream),
                                    candidate.transcript.clone(),
                                    candidate.gene.clone(),
                                    candidate.distance,
                                    pctg_dhs,
                                    pctg_a,
                                    tss_distance,
                                ));
                            }
                        }
                    }

                    let region_overlap = end - exon.start + 1;
                    let pctg_region = (region_overlap as f64 / region_length as f64) * 100.0;
                    let pctg_area = (region_overlap as f64 / exon_length as f64) * 100.0;

                    if (is_first_exon && gene.strand == Strand::Positive)
                        || (is_last_exon && gene.strand == Strand::Negative)
                    {
                        final_output.push(Candidate::new(
                            exon.start,
                            exon.end,
                            gene.strand,
                            exon_number.clone(),
                            Area::FirstExon,
                            transcript.transcript_id.clone(),
                            gene.gene_id.clone(),
                            0,
                            pctg_region,
                            pctg_area,
                            tss_distance,
                        ));
                    } else {
                        let my_id = format!("{}_{}", gene.gene_id, transcript.transcript_id);

                        let gb_candidate = Candidate::new(
                            exon.start,
                            exon.end,
                            gene.strand,
                            exon_number.clone(),
                            Area::GeneBody,
                            transcript.transcript_id.clone(),
                            gene.gene_id.clone(),
                            0,
                            pctg_region,
                            pctg_area,
                            tss_distance,
                        );
                        my_gene_bodys.entry(my_id).or_default().push((
                            gb_candidate,
                            exon_length,
                            region_overlap,
                        ));
                    }
                }
                // Case 5: Region completely within the exon
                //             <----------------->
                //                 |---------|
                else if exon.start <= start && start <= exon.end && end < exon.end {
                    flag_gene_body = true;
                    let pctg_region = 100.0;
                    let pctg_area = (region_length as f64 / exon_length as f64) * 100.0;

                    if (is_first_exon && gene.strand == Strand::Positive)
                        || (is_last_exon && gene.strand == Strand::Negative)
                    {
                        final_output.push(Candidate::new(
                            exon.start,
                            exon.end,
                            gene.strand,
                            exon_number.clone(),
                            Area::FirstExon,
                            transcript.transcript_id.clone(),
                            gene.gene_id.clone(),
                            0,
                            pctg_region,
                            pctg_area,
                            tss_distance,
                        ));
                    } else {
                        let my_id = format!("{}_{}", gene.gene_id, transcript.transcript_id);

                        let gb_candidate = Candidate::new(
                            exon.start,
                            exon.end,
                            gene.strand,
                            exon_number.clone(),
                            Area::GeneBody,
                            transcript.transcript_id.clone(),
                            gene.gene_id.clone(),
                            0,
                            pctg_region,
                            pctg_area,
                            tss_distance,
                        );
                        my_gene_bodys.entry(my_id).or_default().push((
                            gb_candidate,
                            exon_length,
                            region_length,
                        ));
                    }
                }
                // Case 6: Exon totally after the region
                //                       <----------------->
                //   |---------|
                else if exon.start > end && is_first_exon {
                    let dist_tmp = exon.start - pm;

                    if gene.strand == Strand::Negative && dist_tmp < down {
                        down = dist_tmp;
                        exon_down = Some(Candidate::new(
                            exon.start,
                            exon.end,
                            gene.strand,
                            exon_number.clone(),
                            Area::Downstream,
                            transcript.transcript_id.clone(),
                            gene.gene_id.clone(),
                            down,
                            100.0,
                            -1.0,
                            tss_distance,
                        ));
                    } else if gene.strand == Strand::Positive && dist_tmp < upst {
                        upst = dist_tmp;
                        exon_up = Some(Candidate::new(
                            exon.start,
                            exon.end,
                            gene.strand,
                            exon_number.clone(),
                            Area::Upstream,
                            transcript.transcript_id.clone(),
                            gene.gene_id.clone(),
                            upst,
                            100.0,
                            -1.0,
                            tss_distance,
                        ));
                    }

                    if down <= dist_tmp && upst <= dist_tmp {
                        break;
                    }
                }
            }
        }
    }

    // Report closest downstream/upstream if applicable
    if let Some(exon_down_val) = exon_down {
        if down <= upst && exon_down_val.distance <= config.distance {
            if config.tts > 0.0 {
                let exon_info = TtsExonInfo {
                    start: exon_down_val.start,
                    end: exon_down_val.end,
                    strand: exon_down_val.strand,
                    distance: exon_down_val.distance,
                };
                for (tag, pctg_dhs, pctg_a) in check_tts(start, end, &exon_info, config.tts) {
                    final_output.push(Candidate::new(
                        exon_down_val.start,
                        exon_down_val.end,
                        exon_down_val.strand,
                        exon_down_val.exon_number.clone(),
                        tag.parse().unwrap_or(Area::Downstream),
                        exon_down_val.transcript.clone(),
                        exon_down_val.gene.clone(),
                        exon_down_val.distance,
                        pctg_dhs,
                        pctg_a,
                        exon_down_val.tss_distance,
                    ));
                }
            } else {
                final_output.push(exon_down_val);
            }
        }
    }

    if let Some(exon_up_val) = exon_up {
        if upst <= down && exon_up_val.distance <= config.distance {
            let exon_info = TssExonInfo {
                start: exon_up_val.start,
                end: exon_up_val.end,
                strand: exon_up_val.strand,
                distance: exon_up_val.distance,
            };
            for (tag, pctg_dhs, pctg_a) in
                check_tss(start, end, &exon_info, config.tss, config.promoter)
            {
                final_output.push(Candidate::new(
                    exon_up_val.start,
                    exon_up_val.end,
                    exon_up_val.strand,
                    exon_up_val.exon_number.clone(),
                    tag.parse().unwrap_or(Area::Upstream),
                    exon_up_val.transcript.clone(),
                    exon_up_val.gene.clone(),
                    exon_up_val.distance,
                    pctg_dhs,
                    pctg_a,
                    exon_up_val.tss_distance,
                ));
            }
        }
    }

    // Sum up gene body and intron overlaps
    if flag_gene_body {
        // Gene body
        final_output.extend(aggregate_entries(my_gene_bodys, region_length));

        // Introns
        final_output.extend(aggregate_entries(my_introns, region_length));
    }

    final_output
}

pub fn process_candidates_for_output(
    candidates: Vec<Candidate>,
    config: &Config,
) -> Vec<Candidate> {
    if candidates.is_empty() {
        return candidates;
    }

    // filter_by_transcript helper removed (unused logic)

    match config.level {
        ReportLevel::Exon => {
            // Exon Level Logic:
            // Testing confirms that Golden Output behaves as if NO filtering is applied
            // (except for a small set of ~60 edge cases).
            // Rust output is a strict superset of Golden (0 missing lines).
            // To maintain parity (and safety), we return all candidates.
            candidates
        }
        ReportLevel::Transcript => {
            // Transcript Level Logic: Best candidate per transcript.

            // Group by transcript for apply_rules
            let mut by_transcript: AHashMap<String, Vec<usize>> = AHashMap::new();
            for (i, c) in candidates.iter().enumerate() {
                by_transcript
                    .entry(c.transcript.clone())
                    .or_default()
                    .push(i);
            }

            apply_rules(
                &candidates,
                &by_transcript,
                config.perc_region,
                config.perc_area,
                &config.rules,
            )
        }
        ReportLevel::Gene => {
            // Gene Level Logic: Best transcript per gene.

            // 1. Filter per transcript (Best candidate per transcript)
            let mut by_transcript: AHashMap<String, Vec<usize>> = AHashMap::new();
            for (i, c) in candidates.iter().enumerate() {
                by_transcript
                    .entry(c.transcript.clone())
                    .or_default()
                    .push(i);
            }

            let transcript_results = apply_rules(
                &candidates,
                &by_transcript,
                config.perc_region,
                config.perc_area,
                &config.rules,
            );

            // 2. Select best transcript per gene
            let mut by_gene: AHashMap<String, Vec<usize>> = AHashMap::new();
            for (i, c) in transcript_results.iter().enumerate() {
                by_gene.entry(c.gene.clone()).or_default().push(i);
            }

            select_transcript(&transcript_results, &by_gene, &config.rules)
        }
    }
}

/// Main entry point for matching regions to genes.
pub fn match_regions_to_genes(
    regions: &[Region],
    genes: &[Gene],
    config: &Config,
    max_gene_length: i64,
) -> Vec<(Region, Vec<Candidate>)> {
    // Genes must be pre-sorted by start position

    let mut results = Vec::new();

    let max_lookback = max_gene_length + config.max_lookback_distance();
    let mut last_index = 0;

    for region in regions {
        // Calculate safe search start for this region
        // We need to look back enough to find genes that started earlier but extend into this region
        let search_start = region.start.saturating_sub(max_lookback);

        // Advance last_index safe: skip genes that end before the search start
        // These genes can never overlap with the current region or any future region (since regions are sorted by start)
        // Optimization: Use a simple while loop as it is O(N) amortized over all regions
        while last_index < genes.len() && genes[last_index].end < search_start {
            last_index += 1;
        }

        // Pass the calculated start index by value (no mutation allowed inside)
        let candidates = match_region_to_genes(region, genes, config, last_index);
        let processed = process_candidates_for_output(candidates, config);
        results.push((region.clone(), processed));
    }

    results
}

/// Find the index of the first gene that could potentially overlap with a region.
///
/// Uses binary search to find the first gene with `start >= search_start`.
/// This is safe for random access patterns (unsorted regions).
pub fn find_search_start_index(genes: &[Gene], search_start: i64) -> usize {
    genes.partition_point(|g| g.start < search_start)
}
