use rust_htslib::{bam, bam::Read, faidx};
use std::cmp::{max, min};
use std::{collections::HashMap, hash::RandomState};

use log::info;

use crate::args::CommandArgs;
use crate::fasta::index_fasta;
use crate::utils::u8_to_nt;

// We want a lookup table for mapping ascii to nucleotides.
const NTS: [char; 6] = ['A', 'T', 'C', 'G', 'N', '-'];

pub fn bam_parse(args: &CommandArgs) {
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

    println!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        "contig_name",
        "nt",
        "pos",
        "nt_context",
        "depth",
        "ref_cov",
        "rev_cov_abs",
        "deletion_fraction",
        "A",
        "T",
        "C",
        "G",
    );

    // A given position in the reference sequence.
    while let Some(Ok(pileup)) = pileups.next() {
        // Skip low coverage positions.
        let depth = pileup.depth() as usize;

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

        // reference nucleotide.
        let ref_nt = &fasta_info.seq[pos as usize];

        let seq_len = fasta_info.length;

        // For the given position, extract flanking bases to provide context.
        let context_start: usize = max(pos as i32 - args.num_flanking_bases as i32, 0) as usize;
        let context_end: usize =
            min(pos + 1 + (args.num_flanking_bases as u32), seq_len as u32) as usize;
        let context = &fasta_info.seq[context_start..context_end];

        // convert to base.
        let ref_nt_char = u8_to_nt(&ref_nt).unwrap();

        // Hashmap for nt counts.
        let random_state = RandomState::new();
        let mut map: HashMap<char, usize> =
            HashMap::with_capacity_and_hasher(NTS.len(), random_state);

        // Initialize zero count hashmap for nts.
        NTS.iter().for_each(|nt| {
            map.insert(*nt, 0 as usize);
        });

        // Alignment of a single read against a single position.
        for alignment in pileup.alignments() {
            // We just increment num deletions and then skip to next alignment
            if alignment.is_del() || alignment.is_refskip() {
                map.entry('-').and_modify(|c| *c += 1).or_insert(1);
                continue;
            };

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

        // Total reference nucleotide in position.
        let ref_nt_count = map.get(&ref_nt_char).unwrap();

        // Total nucleotides (reference and non-reference).
        let total_nt_count: usize = map.values().sum();

        // Reference coverage only considering nucleotides.
        let ref_frac_cov = *ref_nt_count as f32 / total_nt_count as f32;

        // Reference coverage also considering insertions and deletions.
        let ref_abs_frac_cov = *ref_nt_count as f32 / depth as f32;

        // Calculate deletion fraction (relevant for Nanopore alignments).
        let num_deletions = map.get(&'-').unwrap();
        let deletion_fraction = *num_deletions as f32 / depth as f32;

        if ref_frac_cov < (1.0 - args.min_het_frequency) {
            println!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                reference_name,
                ref_nt_char,
                pos,
                std::str::from_utf8(context).unwrap(),
                depth,
                ref_frac_cov,
                ref_abs_frac_cov,
                deletion_fraction,
                map.get(&'A').unwrap(),
                map.get(&'T').unwrap(),
                map.get(&'C').unwrap(),
                map.get(&'G').unwrap()
            );
        }
    }
}
