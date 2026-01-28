//! CLI entry point for rgmatch.
//!
//! This provides a command-line interface matching the Python implementation.

use anyhow::{bail, Context, Result};
use clap::Parser;
use crossbeam_channel::{bounded, Receiver, Sender};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;
use std::thread;
use std::time::Instant;

use rayon::prelude::*;
use rgmatch::config::Config;
use rgmatch::matcher::overlap::find_search_start_index;
use rgmatch::matcher::{match_region_to_genes, process_candidates_for_output};
use rgmatch::output::{format_output_line, write_header};
use rgmatch::parser::gtf::GtfData;
use rgmatch::parser::{parse_gtf, BedReader};
use rgmatch::types::{Candidate, Region, ReportLevel};

/// Performance metrics for profiling bottlenecks.
/// All times are in nanoseconds.
#[derive(Default)]
struct PerfMetrics {
    /// Total time workers spend on matching (match_region_to_genes + process_candidates)
    worker_matching_ns: AtomicU64,
    /// Total time workers spend waiting to send results on the channel
    worker_channel_wait_ns: AtomicU64,
    /// Total time the writer spends formatting output lines
    writer_format_ns: AtomicU64,
    /// Total time the writer spends on I/O (writeln!)
    writer_io_ns: AtomicU64,
    /// Number of regions processed by workers
    regions_processed: AtomicU64,
    /// Number of output lines written
    lines_written: AtomicU64,
    /// Maximum size of the pending buffer in the writer
    max_pending_size: AtomicU64,
}

impl PerfMetrics {
    fn new() -> Self {
        Self::default()
    }

    fn add_worker_matching(&self, ns: u64) {
        self.worker_matching_ns.fetch_add(ns, Ordering::Relaxed);
    }

    fn add_worker_channel_wait(&self, ns: u64) {
        self.worker_channel_wait_ns.fetch_add(ns, Ordering::Relaxed);
    }

    fn add_writer_format(&self, ns: u64) {
        self.writer_format_ns.fetch_add(ns, Ordering::Relaxed);
    }

    fn add_writer_io(&self, ns: u64) {
        self.writer_io_ns.fetch_add(ns, Ordering::Relaxed);
    }

    fn add_regions_processed(&self, count: u64) {
        self.regions_processed.fetch_add(count, Ordering::Relaxed);
    }

    fn add_lines_written(&self, count: u64) {
        self.lines_written.fetch_add(count, Ordering::Relaxed);
    }

    fn update_max_pending(&self, size: usize) {
        let size = size as u64;
        let mut current = self.max_pending_size.load(Ordering::Relaxed);
        while size > current {
            match self.max_pending_size.compare_exchange_weak(
                current,
                size,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(c) => current = c,
            }
        }
    }

    fn print_summary(&self) {
        let worker_matching_ms = self.worker_matching_ns.load(Ordering::Relaxed) as f64 / 1_000_000.0;
        let worker_channel_wait_ms = self.worker_channel_wait_ns.load(Ordering::Relaxed) as f64 / 1_000_000.0;
        let writer_format_ms = self.writer_format_ns.load(Ordering::Relaxed) as f64 / 1_000_000.0;
        let writer_io_ms = self.writer_io_ns.load(Ordering::Relaxed) as f64 / 1_000_000.0;
        let regions = self.regions_processed.load(Ordering::Relaxed);
        let lines = self.lines_written.load(Ordering::Relaxed);
        let max_pending = self.max_pending_size.load(Ordering::Relaxed);

        eprintln!("\n=== Performance Metrics ===");
        eprintln!("Regions processed: {}", regions);
        eprintln!("Lines written: {}", lines);
        eprintln!();
        eprintln!("Worker time (cumulative across all workers):");
        eprintln!("  Matching:      {:>10.2} ms", worker_matching_ms);
        eprintln!("  Channel wait:  {:>10.2} ms", worker_channel_wait_ms);
        eprintln!();
        eprintln!("Writer time:");
        eprintln!("  Formatting:    {:>10.2} ms", writer_format_ms);
        eprintln!("  I/O:           {:>10.2} ms", writer_io_ms);
        eprintln!();
        eprintln!("Channel congestion:");
        eprintln!("  Max pending results: {} (channel bound: 2000)", max_pending);
        if max_pending >= 1900 {
            eprintln!("  ⚠️  Channel nearly full - WRITER IS BOTTLENECK");
        } else if max_pending < 100 {
            eprintln!("  ✓  Channel uncongested - Workers are bottleneck");
        } else {
            eprintln!("  ~  Moderate congestion - Mixed bottleneck");
        }
        eprintln!();

        // Calculate ratios
        let total_worker = worker_matching_ms + worker_channel_wait_ms;
        let total_writer = writer_format_ms + writer_io_ms;
        if total_worker > 0.0 {
            eprintln!("Worker breakdown:");
            eprintln!("  Matching: {:.1}%", 100.0 * worker_matching_ms / total_worker);
            eprintln!("  Waiting:  {:.1}%", 100.0 * worker_channel_wait_ms / total_worker);
        }
        if total_writer > 0.0 {
            eprintln!("Writer breakdown:");
            eprintln!("  Format: {:.1}%", 100.0 * writer_format_ms / total_writer);
            eprintln!("  I/O:    {:.1}%", 100.0 * writer_io_ms / total_writer);
        }
        eprintln!("=== End Performance Metrics ===\n");
    }
}

/// Genomic region-to-gene matching tool.
///
/// Maps genomic regions from a BED file to gene annotations from a GTF file.
#[derive(Parser, Debug)]
#[command(name = "rgmatch")]
#[command(author, version, about, long_about = None)]
struct Args {
    /// GTF annotation file (required)
    #[arg(short = 'g', long = "gtf")]
    gtf: PathBuf,

    /// Region BED file (required)
    #[arg(short = 'b', long = "bed")]
    bed: PathBuf,

    /// Output file (required)
    #[arg(short = 'o', long = "output")]
    output: PathBuf,

    /// Report level: exon, transcript, or gene
    #[arg(short = 'r', long = "report", default_value = "exon")]
    report: String,

    /// Maximum distance in kb to report associations
    #[arg(short = 'q', long = "distance", default_value = "10")]
    distance: i64,

    /// TSS region distance in bp
    #[arg(short = 't', long = "tss", default_value = "200")]
    tss: i64,

    /// TTS region distance in bp
    #[arg(short = 's', long = "tts", default_value = "0")]
    tts: i64,

    /// Promoter region distance in bp
    #[arg(short = 'p', long = "promoter", default_value = "1300")]
    promoter: i64,

    /// Percentage of the area overlap threshold (0-100)
    #[arg(short = 'v', long = "perc_area", default_value = "90")]
    perc_area: f64,

    /// Percentage of the region overlap threshold (0-100)
    #[arg(short = 'w', long = "perc_region", default_value = "50")]
    perc_region: f64,

    /// Priority rules (comma-separated)
    #[arg(short = 'R', long = "rules", default_value = "TSS,1st_EXON,PROMOTER,TTS,INTRON,GENE_BODY,UPSTREAM,DOWNSTREAM")]
    rules: String,

    /// GTF tag for gene ID
    #[arg(short = 'G', long = "gene", default_value = "gene_id")]
    gene_tag: String,

    /// GTF tag for transcript ID
    #[arg(short = 'T', long = "transcript", default_value = "transcript_id")]
    transcript_tag: String,

    /// Number of worker threads (0 = auto-detect, 1 = sequential)
    #[arg(long = "threads", short = 'j', default_value = "8")]
    threads: usize,

    /// Batch size for streaming BED regions
    #[arg(long = "batch-size", default_value = "5000")]
    batch_size: usize,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Validate inputs
    if !args.gtf.exists() {
        bail!("GTF file not found: {}", args.gtf.display());
    }
    if !args.bed.exists() {
        bail!("BED file not found: {}", args.bed.display());
    }

    // Parse report level
    let level = ReportLevel::from_str(&args.report).context(
        "Report can only be one of the following: exon, transcript or gene",
    )?;

    // Build configuration
    let mut config = Config::new();
    config.level = level;

    // Set distance (convert from kb to bp)
    if args.distance >= 0 {
        config.set_distance_kb(args.distance);
    }

    // Set TSS distance
    if args.tss >= 0 {
        config.tss = args.tss as f64;
    } else {
        bail!("The TSS distance cannot be lower than 0 bps.");
    }

    // Set TTS distance
    if args.tts >= 0 {
        config.tts = args.tts as f64;
    } else {
        bail!("The TTS distance cannot be lower than 0 bps.");
    }

    // Set promoter distance
    if args.promoter >= 0 {
        config.promoter = args.promoter as f64;
    } else {
        bail!("The promoter distance cannot be lower than 0 bps.");
    }

    // Set percentage thresholds
    if args.perc_area >= 0.0 && args.perc_area <= 100.0 {
        config.perc_area = args.perc_area;
    } else {
        bail!("The percentage of area defined was wrong. It should range between 0 and 100.");
    }

    if args.perc_region >= 0.0 && args.perc_region <= 100.0 {
        config.perc_region = args.perc_region;
    } else {
        bail!("The percentage of region defined was wrong. It should range between 0 and 100.");
    }

    // Parse rules
    if !config.parse_rules(&args.rules) {
        bail!("Rules not properly passed.");
    }

    // Set GTF tags
    config.gene_id_tag = args.gene_tag.clone();
    config.transcript_id_tag = args.transcript_tag.clone();

    // Parse GTF file
    eprintln!("Parsing GTF file: {}", args.gtf.display());
    let mut gtf_data = parse_gtf(&args.gtf, &config.gene_id_tag, &config.transcript_id_tag)?;

    // Pre-sort genes for deterministic matching and performance
    gtf_data.genes_by_chrom.values_mut().collect::<Vec<_>>().par_iter_mut().for_each(|genes| {
        genes.sort_by(|a, b| a.start.cmp(&b.start).then(a.gene_id.cmp(&b.gene_id)));
    });

    // Validate batch_size
    if args.batch_size == 0 {
        bail!("Batch size must be greater than 0");
    }

    // Determine thread count
    let num_threads = if args.threads == 0 {
        num_cpus::get()
    } else {
        args.threads
    };

    if num_threads == 1 {
        // Use original sequential implementation
        run_sequential(&args, &gtf_data, &config)?;
    } else {
        // Use parallel pipeline
        run_parallel(&args, gtf_data, &config, num_threads)?;
    }

    eprintln!("Done!");
    Ok(())
}

/// Sequential implementation with streaming.
fn run_sequential(args: &Args, gtf_data: &GtfData, config: &Config) -> Result<()> {
    eprintln!("Processing BED file: {}", args.bed.display());
    
    // Initialize streaming reader
    let mut bed_reader = BedReader::new(&args.bed)?;
    
    // Output writer
    eprintln!("Writing output to: {}", args.output.display());
    let file = File::create(&args.output).context("Failed to create output file")?;
    let mut writer = BufWriter::new(file);

    let mut header_written = false;
    
    // Optimization state
    let mut last_chrom = String::new();
    let mut last_start = -1;
    let mut last_index = 0;

    // Process in chunks
    while let Some(chunk) = bed_reader.read_chunk(args.batch_size)? {
        if !header_written {
            let num_meta = bed_reader.num_meta_columns();
            write_header(&mut writer, num_meta)?;
            header_written = true;
        }

        for region in chunk {
            // Find genes for chrom
            if let Some(genes) = gtf_data.genes_by_chrom.get(&region.chrom) {
                 let max_len = *gtf_data.max_lengths.get(&region.chrom).unwrap_or(&0);
                 
                 // Calculate safe search start (region start - max_len - distance)
                 // Note: we must match the logic in match_regions_to_genes regarding max_lookback
                 let max_lookback = max_len + config.max_lookback_distance();
                 let search_start = region.start.saturating_sub(max_lookback);
                 
                 let start_index;
                 if region.chrom == last_chrom && region.start >= last_start {
                     // Optimistic: advance from last_index
                     let mut idx = last_index;
                     // Skip genes that end before search_start
                     while idx < genes.len() && genes[idx].end < search_start {
                         idx += 1;
                     }
                     start_index = idx;
                 } else {
                     // Reset / Binary search
                     start_index = find_search_start_index(genes, search_start);
                 }
                 
                 // Update cache
                 last_chrom = region.chrom.clone();
                 last_start = region.start;
                 last_index = start_index;
                 
                 // Match
                 let candidates = match_region_to_genes(&region, genes, config, start_index);
                 let processed = process_candidates_for_output(candidates, config);
                 
                 // Write line
                 for candidate in processed {
                     let line = format_output_line(&region, &candidate);
                     writeln!(writer, "{}", line)?;
                 }
            } else {
                // If chromosome not in GTF, verify if we should reset cache?
                // Probably yes to be safe, though chrom changed so next valid chrom will trigger binary search.
                last_chrom = region.chrom.clone();
            }
        }
    }
    
    if !header_written {
         // File was empty
         write_header(&mut writer, 0)?;
    }
    
    writer.flush()?;
    Ok(())
}

/// Work item for the parallel pipeline.
struct WorkItem {
    /// Sequence number for ordering (file order).
    seq_id: u64,
    /// Regions to process (all from same chromosome, in file order).
    regions: Vec<Region>,
}

/// Result from processing a work item.
struct WorkResult {
    /// Sequence number matching the input WorkItem.
    seq_id: u64,
    /// Processing results in the same order as input regions.
    results: Vec<(Region, Vec<Candidate>)>,
}

/// Parallel implementation using per-chromosome work distribution.
///
/// To ensure byte-for-byte compatibility with sequential mode, we:
/// 1. Parse the entire BED file and group regions by chromosome
/// 2. Distribute chromosomes to workers (each chromosome is one work item)
/// 3. Write results in sorted chromosome order
/// Parallel implementation with streaming.
fn run_parallel(
    args: &Args,
    gtf_data: GtfData,
    config: &Config,
    num_threads: usize,
) -> Result<()> {
    eprintln!("Using parallel mode with {} threads", num_threads);

    // Create performance metrics
    let metrics = Arc::new(PerfMetrics::new());

    // Create channels
    let (work_tx, work_rx): (Sender<WorkItem>, Receiver<WorkItem>) = bounded(100);
    // Increased buffer for results to avoid blocking workers
    let (result_tx, result_rx): (Sender<WorkResult>, Receiver<WorkResult>) = bounded(2000);

    // Shared GTF data for workers
    let gtf_arc = Arc::new(gtf_data);
    let config_arc = Arc::new(config.clone());

    // Spawn writer thread
    let output_path = args.output.clone();

    let (header_tx, header_rx) = bounded(1);

    let writer_handle = thread::spawn({
        let result_rx = result_rx.clone();
        let metrics = Arc::clone(&metrics);
        move || -> Result<usize> { write_results_ordered(&output_path, result_rx, header_rx, &metrics) }
    });

    // Spawn worker threads using rayon's thread pool
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .context("Failed to create thread pool")?;

    // Clone references for the worker scope
    let gtf_for_workers = Arc::clone(&gtf_arc);
    let config_for_workers = Arc::clone(&config_arc);
    let work_rx_for_workers = work_rx.clone();
    let result_tx_for_workers = result_tx.clone();
    let metrics_for_workers = Arc::clone(&metrics);

    // Spawn workers in a separate thread to avoid blocking
    let workers_handle = thread::spawn(move || {
        pool.scope(|s| {
            for _ in 0..num_threads {
                let work_rx = work_rx_for_workers.clone();
                let result_tx = result_tx_for_workers.clone();
                let gtf = Arc::clone(&gtf_for_workers);
                let cfg = Arc::clone(&config_for_workers);
                let metrics = Arc::clone(&metrics_for_workers);

                s.spawn(move |_| {
                    worker_loop(work_rx, result_tx, gtf, cfg, &metrics);
                });
            }
        });
    });

    // Producer: Read BED in chunks
    eprintln!("Processing BED file: {}", args.bed.display());
    let mut bed_reader = BedReader::new(&args.bed)?;
    
    let mut global_seq_id = 0;
    
    // Send header info immediately if possible? No, header depends on first line read usually.
    // BedReader logic: read_chunk updates num_meta_columns.
    // So we need to read first chunk.
    
    loop {
        match bed_reader.read_chunk(args.batch_size)? {
            Some(chunk) => {
                if global_seq_id == 0 {
                    // Send header info
                    let _ = header_tx.send(bed_reader.num_meta_columns());
                }
                
                let work_item = WorkItem {
                    seq_id: global_seq_id,
                    regions: chunk,
                };
                
                if work_tx.send(work_item).is_err() {
                    break;
                }
                global_seq_id += 1;
            },
            None => break,
        }
    }
    
    // If loop finished and global_seq_id is 0, file was empty.
    if global_seq_id == 0 {
        let _ = header_tx.send(0);
    }

    // Close work channel to signal workers to exit
    drop(work_tx);
    drop(header_tx); // Close header channel too

    // Wait for workers to finish
    workers_handle
        .join()
        .map_err(|_| anyhow::anyhow!("Worker thread panicked"))?;

    // Close result channel to signal writer to finish
    drop(result_tx);

    // Wait for writer and get the results
    let lines_written = writer_handle
        .join()
        .map_err(|_| anyhow::anyhow!("Writer thread panicked"))??;

    eprintln!(
        "Writing output to: {} ({} lines)",
        args.output.display(),
        lines_written
    );

    // Print performance metrics
    metrics.print_summary();

    Ok(())
}

/// Worker loop: receives work items and sends results.
fn worker_loop(
    work_rx: Receiver<WorkItem>,
    result_tx: Sender<WorkResult>,
    gtf: Arc<GtfData>,
    config: Arc<Config>,
    metrics: &PerfMetrics,
) {
    // Optimization state per worker
    let mut last_chrom = String::new();
    let mut last_start = -1;
    let mut last_index = 0;

    while let Ok(work_item) = work_rx.recv() {
        let num_regions = work_item.regions.len() as u64;

        // Time the matching work
        let match_start = Instant::now();
        let results = process_work_item(&work_item, &gtf, &config, &mut last_chrom, &mut last_start, &mut last_index);
        let match_elapsed = match_start.elapsed();
        metrics.add_worker_matching(match_elapsed.as_nanos() as u64);
        metrics.add_regions_processed(num_regions);

        let work_result = WorkResult {
            seq_id: work_item.seq_id,
            results,
        };

        // Time the channel send (how long we wait if channel is full)
        let send_start = Instant::now();
        let send_result = result_tx.send(work_result);
        let send_elapsed = send_start.elapsed();
        metrics.add_worker_channel_wait(send_elapsed.as_nanos() as u64);

        if send_result.is_err() {
            break;
        }
    }
}

/// Process a single work item (a chunk of regions).
fn process_work_item(
    work_item: &WorkItem,
    gtf: &GtfData,
    config: &Config,
    last_chrom: &mut String,
    last_start: &mut i64,
    last_index: &mut usize,
) -> Vec<(Region, Vec<Candidate>)> {
    let mut results = Vec::with_capacity(work_item.regions.len());

    for region in &work_item.regions {
        if let Some(genes) = gtf.genes_by_chrom.get(&region.chrom) {
             let max_len = *gtf.max_lengths.get(&region.chrom).unwrap_or(&0);
             
             let max_lookback = max_len + config.max_lookback_distance();
             let search_start = region.start.saturating_sub(max_lookback);
             
             let start_index;
             if *last_chrom == region.chrom && region.start >= *last_start {
                 let mut idx = *last_index;
                 while idx < genes.len() && genes[idx].end < search_start {
                     idx += 1;
                 }
                 start_index = idx;
             } else {
                 start_index = find_search_start_index(genes, search_start);
             }
             
             *last_chrom = region.chrom.clone();
             *last_start = region.start;
             *last_index = start_index;
             
             let candidates = match_region_to_genes(region, genes, config, start_index);
             let processed = process_candidates_for_output(candidates, config);
             results.push((region.clone(), processed));
        } else {
             // Chromosome not found, but we must record it in output as processed (with empty candidates) 
             // wait, match_region_to_genes returns Vec<Candidate>.
             // If no genes, results is empty.
             // But original code: "If the current gene also covers the region..."
             // If no genes for chrom, we skip?
             // run_sequential (original) printd "Warning: {} not found" and skipped.
             // We should maintain parity.
             // If skipping, we don't push to results?
             // But we need to maintain order?
             // Actually, if a region has no matches, it produces no output lines. 
             // So skipping here is fine.
             *last_chrom = region.chrom.clone();
        }
    }
    
    results
}

/// Write results in order, buffering out-of-order results.
fn write_results_ordered(
    output_path: &PathBuf,
    result_rx: Receiver<WorkResult>,
    header_rx: Receiver<usize>,
    metrics: &PerfMetrics,
) -> Result<usize> {
    let file = File::create(output_path).context("Failed to create output file")?;
    let mut writer = BufWriter::new(file);

    // Get header info (blocking until first chunk read or empty file)
    let num_meta_columns = header_rx.recv().unwrap_or(0);
    write_header(&mut writer, num_meta_columns)?;

    // Buffer for out-of-order results
    let mut pending: BTreeMap<u64, WorkResult> = BTreeMap::new();
    let mut next_expected: u64 = 0;
    let mut lines_written: usize = 0;

    for result in result_rx {
        pending.insert(result.seq_id, result);

        // Track max pending size for congestion analysis
        metrics.update_max_pending(pending.len());

        // Write all ready consecutive results
        while let Some(r) = pending.remove(&next_expected) {
            for (region, candidates) in &r.results {
                for candidate in candidates {
                    // Time formatting
                    let format_start = Instant::now();
                    let line = format_output_line(region, candidate);
                    let format_elapsed = format_start.elapsed();
                    metrics.add_writer_format(format_elapsed.as_nanos() as u64);

                    // Time I/O
                    let io_start = Instant::now();
                    writeln!(writer, "{}", line)?;
                    let io_elapsed = io_start.elapsed();
                    metrics.add_writer_io(io_elapsed.as_nanos() as u64);

                    lines_written += 1;
                }
            }
            next_expected += 1;
        }
    }

    metrics.add_lines_written(lines_written as u64);
    writer.flush()?;
    Ok(lines_written)
}
