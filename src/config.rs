//! Configuration and defaults for rgmatch.
//!
//! This module contains the configuration structure and default values
//! that control the region-to-gene matching behavior.

use crate::types::{Area, ReportLevel};

/// Default rules priority order.
pub const DEFAULT_RULES: [Area; 8] = [
    Area::Tss,
    Area::FirstExon,
    Area::Promoter,
    Area::Tts,
    Area::Intron,
    Area::GeneBody,
    Area::Upstream,
    Area::Downstream,
];

/// Configuration for the region-to-gene matching process.
#[derive(Debug, Clone)]
pub struct Config {
    /// Priority rules for resolving ties.
    pub rules: Vec<Area>,
    /// Percentage of the area overlapped threshold.
    pub perc_area: f64,
    /// Percentage of the region overlapped threshold.
    pub perc_region: f64,
    /// TSS region distance in bp.
    pub tss: f64,
    /// TTS region distance in bp.
    pub tts: f64,
    /// Promoter region distance in bp.
    pub promoter: f64,
    /// Maximum distance to report associations in bp.
    pub distance: i64,
    /// Report level (exon, transcript, or gene).
    pub level: ReportLevel,
    /// GTF tag for gene ID.
    pub gene_id_tag: String,
    /// GTF tag for transcript ID.
    pub transcript_id_tag: String,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            rules: DEFAULT_RULES.to_vec(),
            perc_area: 90.0,
            perc_region: 50.0,
            tss: 200.0,
            tts: 0.0,
            promoter: 1300.0,
            distance: 10000, // 10kb default (stored in bp)
            level: ReportLevel::Exon,
            gene_id_tag: "gene_id".to_string(),
            transcript_id_tag: "transcript_id".to_string(),
        }
    }
}

impl Config {
    /// Create a new config with default values.
    pub fn new() -> Self {
        Self::default()
    }

    /// Parse and validate priority rules from a comma-separated string.
    ///
    /// Returns true if all 8 valid tags were provided, false otherwise.
    pub fn parse_rules(&mut self, rules_str: &str) -> bool {
        let valid_tags = [
            "TSS",
            "1st_EXON",
            "PROMOTER",
            "TTS",
            "INTRON",
            "GENE_BODY",
            "UPSTREAM",
            "DOWNSTREAM",
        ];

        let mut new_rules = Vec::new();
        let parts: Vec<&str> = rules_str.split(',').collect();

        for tag in parts {
            if valid_tags.contains(&tag) {
                if let Some(area) = Area::from_str(tag) {
                    // Only add if not already present
                    if !new_rules.contains(&area) {
                        new_rules.push(area);
                    }
                }
            }
        }

        if new_rules.len() == 8 {
            self.rules = new_rules;
            true
        } else {
            false
        }
    }

    /// Set distance in kb (converts to bp internally).
    pub fn set_distance_kb(&mut self, kb: i64) {
        if kb >= 0 {
            self.distance = kb * 1000;
        }
    }

    /// Get the maximum distance to consider for lookback
    pub fn max_lookback_distance(&self) -> i64 {
        let max_float = self.tss.max(self.tts).max(self.promoter);
        self.distance.max(max_float as i64)
    }
}

#[cfg(test)]
mod tests {
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
        assert_eq!(config.level, ReportLevel::Exon);
        assert_eq!(config.gene_id_tag, "gene_id");
        assert_eq!(config.transcript_id_tag, "transcript_id");
    }

    #[test]
    fn test_parse_rules_valid() {
        let mut config = Config::new();
        let result = config.parse_rules(
            "DOWNSTREAM,UPSTREAM,GENE_BODY,INTRON,TTS,PROMOTER,1st_EXON,TSS",
        );
        assert!(result);
        assert_eq!(config.rules.len(), 8);
        assert_eq!(config.rules[0], Area::Downstream);
        assert_eq!(config.rules[7], Area::Tss);
    }

    #[test]
    fn test_parse_rules_default_order() {
        let mut config = Config::new();
        let result = config.parse_rules(
            "TSS,1st_EXON,PROMOTER,TTS,INTRON,GENE_BODY,UPSTREAM,DOWNSTREAM",
        );
        assert!(result);
        assert_eq!(config.rules.len(), 8);
    }

    #[test]
    fn test_parse_rules_missing_tags() {
        let mut config = Config::new();
        let result = config.parse_rules("TSS,1st_EXON,PROMOTER");
        assert!(!result);
    }

    #[test]
    fn test_parse_rules_unknown_tag() {
        let mut config = Config::new();
        let result = config.parse_rules(
            "TSS,1st_EXON,PROMOTER,TTS,INTRON,GENE_BODY,UPSTREAM,UNKNOWN",
        );
        assert!(!result);
    }

    #[test]
    fn test_parse_rules_duplicates() {
        let mut config = Config::new();
        let result = config.parse_rules("TSS,TSS,TSS,TSS,TSS,TSS,TSS,TSS");
        assert!(!result);
    }

    #[test]
    fn test_parse_rules_case_sensitive() {
        let mut config = Config::new();
        let result = config.parse_rules(
            "tss,1st_exon,promoter,tts,intron,gene_body,upstream,downstream",
        );
        assert!(!result);
    }

    #[test]
    fn test_parse_rules_empty() {
        let mut config = Config::new();
        let result = config.parse_rules("");
        assert!(!result);
    }

    #[test]
    fn test_parse_rules_whitespace() {
        let mut config = Config::new();
        let result = config.parse_rules(
            "TSS, 1st_EXON, PROMOTER, TTS, INTRON, GENE_BODY, UPSTREAM, DOWNSTREAM",
        );
        assert!(!result); // Spaces make tags invalid
    }

    #[test]
    fn test_set_distance_kb() {
        let mut config = Config::new();
        config.set_distance_kb(20);
        assert_eq!(config.distance, 20000);

        config.set_distance_kb(-1);
        assert_eq!(config.distance, 20000); // Should not change for negative values
    }
}
