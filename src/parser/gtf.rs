//! GTF file parser with gzip support.
//!
//! Parses GTF (Gene Transfer Format) annotation files to build a hierarchical
//! structure of genes, transcripts, and exons organized by chromosome.

use ahash::AHashMap;
use anyhow::{Context, Result};
use std::fs::File;
use std::io::BufRead;
use std::path::Path;

use crate::parser::util::create_buffered_reader;
use crate::types::{Exon, Gene, Strand, Transcript};

/// Result of parsing a GTF file.
#[derive(Clone)]
pub struct GtfData {
    /// Genes organized by chromosome.
    pub genes_by_chrom: AHashMap<String, Vec<Gene>>,
    /// Maximum gene length per chromosome.
    pub max_lengths: AHashMap<String, i64>,
}

/// Parse a GTF file and return organized gene data.
///
/// Supports both plain text and gzip-compressed GTF files.
pub fn parse_gtf(path: &Path, gene_id_tag: &str, transcript_id_tag: &str) -> Result<GtfData> {
    let file = File::open(path).context("Failed to open GTF file")?;
    let reader = create_buffered_reader(file, path);

    parse_gtf_reader(reader, gene_id_tag, transcript_id_tag)
}

/// Parse GTF data from a reader.
fn parse_gtf_reader<R: BufRead>(
    reader: R,
    gene_id_tag: &str,
    transcript_id_tag: &str,
) -> Result<GtfData> {
    // Maps to track all genes and transcripts
    let mut all_genes: AHashMap<String, Gene> = AHashMap::new();
    let mut all_transcripts: AHashMap<String, usize> = AHashMap::new(); // transcript_id -> index in gene
    let mut gene_to_transcripts: AHashMap<String, Vec<String>> = AHashMap::new(); // gene_id -> transcript_ids

    // Genes organized by chromosome
    let mut genes_by_chrom: AHashMap<String, Vec<String>> = AHashMap::new(); // chrom -> gene_ids (in order added)

    // Flags to track if transcript and gene entries exist in GTF
    let mut gene_flag = false;
    let mut trans_flag = false;

    for line_result in reader.lines() {
        let line = line_result.context("Failed to read GTF line")?;

        // Skip empty lines and comments
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let chrom = fields[0];
        let feature_type = fields[2];
        let start: i64 = fields[3]
            .parse()
            .context("Failed to parse start coordinate")?;
        let end: i64 = fields[4]
            .parse()
            .context("Failed to parse end coordinate")?;
        let strand_str = fields[6];
        let attributes = fields[8];

        let strand = match strand_str.parse::<Strand>() {
            Ok(s) => s,
            Err(_) => continue, // Skip entries without valid strand
        };

        match feature_type {
            "exon" => {
                let gene_id = extract_attribute(attributes, gene_id_tag)
                    .context("Failed to extract gene_id from exon")?;
                let transcript_id = extract_attribute(attributes, transcript_id_tag)
                    .context("Failed to extract transcript_id from exon")?;

                // Create or get gene
                if !all_genes.contains_key(&gene_id) {
                    all_genes.insert(gene_id.clone(), Gene::new(gene_id.clone(), strand));
                    genes_by_chrom
                        .entry(chrom.to_string())
                        .or_default()
                        .push(gene_id.clone());
                }

                // Create or get transcript
                let is_new_transcript = !all_transcripts.contains_key(&transcript_id);
                if is_new_transcript {
                    let gene = all_genes.get_mut(&gene_id).unwrap();
                    let transcript_idx = gene.transcripts.len();
                    gene.add_transcript(Transcript::new(transcript_id.clone()));
                    all_transcripts.insert(transcript_id.clone(), transcript_idx);
                    gene_to_transcripts
                        .entry(gene_id.clone())
                        .or_default()
                        .push(transcript_id.clone());
                }

                // Add exon to transcript
                let exon = Exon::new(start, end);
                let transcript_idx = all_transcripts[&transcript_id];
                let gene = all_genes.get_mut(&gene_id).unwrap();
                gene.transcripts[transcript_idx].add_exon(exon);
            }
            "transcript" => {
                trans_flag = true;

                let gene_id = extract_attribute(attributes, gene_id_tag)
                    .context("Failed to extract gene_id from transcript")?;
                let transcript_id = extract_attribute(attributes, transcript_id_tag)
                    .context("Failed to extract transcript_id from transcript")?;

                // Create or get gene
                if !all_genes.contains_key(&gene_id) {
                    all_genes.insert(gene_id.clone(), Gene::new(gene_id.clone(), strand));
                    genes_by_chrom
                        .entry(chrom.to_string())
                        .or_default()
                        .push(gene_id.clone());
                }

                // Create or get transcript
                let is_new_transcript = !all_transcripts.contains_key(&transcript_id);
                if is_new_transcript {
                    let gene = all_genes.get_mut(&gene_id).unwrap();
                    let transcript_idx = gene.transcripts.len();
                    gene.add_transcript(Transcript::new(transcript_id.clone()));
                    all_transcripts.insert(transcript_id.clone(), transcript_idx);
                    gene_to_transcripts
                        .entry(gene_id.clone())
                        .or_default()
                        .push(transcript_id.clone());
                }

                // Set transcript boundaries
                let transcript_idx = all_transcripts[&transcript_id];
                let gene = all_genes.get_mut(&gene_id).unwrap();
                gene.transcripts[transcript_idx].set_length(start, end);
            }
            "gene" => {
                gene_flag = true;

                let gene_id = extract_attribute(attributes, gene_id_tag)
                    .context("Failed to extract gene_id from gene")?;

                // Create or get gene
                if !all_genes.contains_key(&gene_id) {
                    all_genes.insert(gene_id.clone(), Gene::new(gene_id.clone(), strand));
                    genes_by_chrom
                        .entry(chrom.to_string())
                        .or_default()
                        .push(gene_id.clone());
                }

                // Set gene boundaries
                all_genes.get_mut(&gene_id).unwrap().set_length(start, end);
            }
            _ => {
                // Skip other feature types
            }
        }
    }

    // Post-processing: check exon numbers and calculate sizes
    for gene in all_genes.values_mut() {
        let strand = gene.strand;
        for transcript in &mut gene.transcripts {
            // Renumber exons based on strand
            transcript.renumber_exons(strand);

            // Calculate transcript size if not set from transcript entry
            if !trans_flag {
                transcript.calculate_size();
            }
        }
    }

    // Calculate gene sizes if not set from gene entries
    if !gene_flag {
        for gene in all_genes.values_mut() {
            gene.calculate_size();
        }
    }

    // Build final genes_by_chrom with actual Gene objects
    let mut result_genes: AHashMap<String, Vec<Gene>> = AHashMap::new();
    let mut max_lengths: AHashMap<String, i64> = AHashMap::new();

    for (chrom, gene_ids) in genes_by_chrom {
        let genes: Vec<Gene> = gene_ids
            .into_iter()
            .filter_map(|id| all_genes.remove(&id))
            .collect();

        let max_len = genes.iter().map(|g| g.end - g.start).max().unwrap_or(0);
        max_lengths.insert(chrom.clone(), max_len);

        result_genes.insert(chrom, genes);
    }

    Ok(GtfData {
        genes_by_chrom: result_genes,
        max_lengths,
    })
}

/// Extract an attribute value from the GTF attributes string.
///
/// GTF attributes are in the format: key "value"; key "value"; ...
fn extract_attribute(attributes: &str, key: &str) -> Option<String> {
    // Find the key
    let key_pattern = format!("{} ", key);
    let start_idx = attributes.find(&key_pattern)?;

    // Find the value between quotes
    let after_key = &attributes[start_idx + key_pattern.len()..];

    // Find first quote
    let first_quote = after_key.find('"')?;
    let after_first_quote = &after_key[first_quote + 1..];

    // Find second quote
    let second_quote = after_first_quote.find('"')?;

    Some(after_first_quote[..second_quote].to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    #[test]
    fn test_extract_attribute() {
        let attrs = r#"gene_id "ENSG00000279493.1"; transcript_id "ENST00000624081.1"; gene_type "artifact";"#;

        assert_eq!(
            extract_attribute(attrs, "gene_id"),
            Some("ENSG00000279493.1".to_string())
        );
        assert_eq!(
            extract_attribute(attrs, "transcript_id"),
            Some("ENST00000624081.1".to_string())
        );
        assert_eq!(
            extract_attribute(attrs, "gene_type"),
            Some("artifact".to_string())
        );
        assert_eq!(extract_attribute(attrs, "nonexistent"), None);
    }

    #[test]
    fn test_parse_gtf_reader() {
        let gtf_content = r#"##description: test
chr1	TEST	gene	1000	2000	.	+	.	gene_id "G1"; gene_name "Gene1";
chr1	TEST	transcript	1000	2000	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	TEST	exon	1000	1200	.	+	.	gene_id "G1"; transcript_id "T1"; exon_number 1;
chr1	TEST	exon	1500	2000	.	+	.	gene_id "G1"; transcript_id "T1"; exon_number 2;
"#;

        let reader = BufReader::new(gtf_content.as_bytes());
        let result = parse_gtf_reader(reader, "gene_id", "transcript_id").unwrap();

        assert!(result.genes_by_chrom.contains_key("chr1"));
        let genes = &result.genes_by_chrom["chr1"];
        assert_eq!(genes.len(), 1);

        let gene = &genes[0];
        assert_eq!(gene.gene_id, "G1");
        assert_eq!(gene.strand, Strand::Positive);
        assert_eq!(gene.transcripts.len(), 1);

        let transcript = &gene.transcripts[0];
        assert_eq!(transcript.transcript_id, "T1");
        assert_eq!(transcript.exons.len(), 2);

        // Check exon numbering for positive strand
        assert_eq!(transcript.exons[0].start, 1000);
        assert_eq!(transcript.exons[0].exon_number, Some("1".to_string()));
        assert_eq!(transcript.exons[1].start, 1500);
        assert_eq!(transcript.exons[1].exon_number, Some("2".to_string()));
    }

    #[test]
    fn test_parse_gtf_negative_strand() {
        let gtf_content = r#"chr1	TEST	exon	1000	1200	.	-	.	gene_id "G1"; transcript_id "T1";
chr1	TEST	exon	1500	2000	.	-	.	gene_id "G1"; transcript_id "T1";
"#;

        let reader = BufReader::new(gtf_content.as_bytes());
        let result = parse_gtf_reader(reader, "gene_id", "transcript_id").unwrap();

        let gene = &result.genes_by_chrom["chr1"][0];
        let transcript = &gene.transcripts[0];

        // For negative strand: first (lowest) gets N, last (highest) gets 1
        assert_eq!(transcript.exons[0].start, 1000);
        assert_eq!(transcript.exons[0].exon_number, Some("2".to_string()));
        assert_eq!(transcript.exons[1].start, 1500);
        assert_eq!(transcript.exons[1].exon_number, Some("1".to_string()));
    }
}
