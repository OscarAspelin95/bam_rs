use rust_htslib::{bam, bam::Read, faidx};
use std::cmp::{max, min};
use std::io::Write;
use std::{collections::HashMap, hash::RandomState};

use log::info;
use simple_logger::SimpleLogger;

mod utils;
use utils::{find_homopolymers, u8_to_nt};

mod args;
use args::{get_args, CommandArgs};

mod fasta;
use fasta::index_fasta;

// We want a lookup table for mapping ascii to nucleotides.
static NTS: [char; 6] = ['A', 'T', 'C', 'G', 'N', '-'];

fn bam_parse(args: &CommandArgs) {
    // Bam reading
    let mut bam = bam::Reader::from_path(&args.bam).expect("Failed to read BAM file {bam_test}.");
    let mut pileups = bam.pileup();

    // We need to read the indexed fasta
    info!("Parsing indexed fasta file...");
    let faidx_reader =
        faidx::Reader::from_path(&args.fasta).expect("Failed to read fasta file {fasta_test}.");

    let seq_names: Vec<String> = faidx_reader
        .seq_names()
        .expect("Failed to extract reference names from {fasta_test}.");

    let fasta_index = index_fasta(&faidx_reader, &seq_names);

    // Consider using next_chunks here if possible (probably not), will probably
    // enable running with Rayon for parallel processing.

    // A given position in the reference sequence.
    while let Some(Ok(pileup)) = pileups.next() {
        let depth = pileup.depth() as usize;

        // Skip low coverage regions.
        if depth < args.min_read_depth as usize {
            continue;
        }

        // Reference related.
        let pos = pileup.pos();
        let reference_id = pileup.tid();

        // Extract reference nucleotide at this position.
        // NOTE - It is yet unknown what is the most efficient way to extract a nt from a given
        // reference name and position. The current implementation has a pre-constructed hashmap
        // of structure {"name"; seq}, for which we extract the nt per reference name and position.
        let reference_name = faidx_reader
            .seq_name(reference_id as i32)
            .expect("Failed to extract reference name from reference id {reference_id}");

        let fasta_info = fasta_index.get(reference_name.as_str()).unwrap();

        // reference nucleotide as u8
        // NOTE - possible improvement - extract context window and let ref_nt = middle value.
        let ref_nt = &fasta_info.seq[pos as usize];

        //
        let seq_len = fasta_info.length;

        let context_start: usize = max(pos as i32 - args.num_flanking_bases as i32, 0) as usize;
        let context_end = min(pos + 1 + (args.num_flanking_bases as u32), seq_len as u32) as usize;

        let context = &fasta_info.seq[context_start..context_end];

        // convert to base.
        let s = u8_to_nt(&ref_nt).unwrap();

        // Hashmap for nt counts.
        let random_state = RandomState::new();
        let mut map = HashMap::with_capacity_and_hasher(NTS.len(), random_state);

        // Initialize zero count hashmap for nts.
        NTS.iter().for_each(|nt| {
            map.insert(nt, 0 as usize);
        });

        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            "contig_name", "nt", "pos", "nt_context", "depth", "ref_cov", "deletion_fraction",
        );

        // Alignment of a single read against a single position.
        for alignment in pileup.alignments() {
            // We just increment num deletions and then skip to next alignment
            if alignment.is_del() || alignment.is_refskip() {
                match map.get_mut(&'-') {
                    Some(x) => *x += 1,
                    _ => panic!("Deletion key does not exist in lookup table."),
                }
                continue;
            }

            // Position in read.
            let read_pos = alignment
                .qpos()
                .expect("Failed to extract read position in read position.");

            // read base
            let base = alignment.record().seq()[read_pos];

            let read_nt = u8_to_nt(&base).unwrap();

            match map.get_mut(&read_nt) {
                Some(x) => *x += 1,
                _ => panic!("Invalid nt lookup: {read_nt}"),
            }
        }

        // Now, we can also extract the nt fraction. E.g., if the base is A,
        // what percentage of reads had an A in this position?
        // NOTE - here we can calculate it in two different ways:
        // * Either based on read depth (will include deletions).
        // * Or, based on e.g., count of A / sum(count of A/T/C/G/N) (exclude deletions).

        let ref_nt_count = map.get(&s).unwrap();

        let total_nt_count: usize = map.values().sum();

        // Only considering nucleotides A/T/C/G/N
        let ref_frac_cov = *ref_nt_count as f32 / total_nt_count as f32;

        // Considering total read depth (including deletions, etc)
        // let ref_abs_frac_cov = *ref_nt_count as f32 / depth as f32;

        let num_deletions = map.get(&'-').unwrap();

        let deletion_fraction = *num_deletions as f32 / depth as f32;

        if ref_frac_cov < args.min_het_frequency {
            println!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                reference_name,
                s,
                pos,
                String::from_utf8(context.to_vec()).unwrap(),
                depth,
                ref_frac_cov,
                deletion_fraction,
            );
            std::io::stdout().flush().unwrap();
        }
    }
}

fn main() {
    SimpleLogger::new().init().unwrap();
    let args = get_args();

    bam_parse(&args);
}
