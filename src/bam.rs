use log::info;
use rust_htslib::{bam, bam::Read, faidx};
use std::cmp::{max, min};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::{collections::HashMap, hash::RandomState};

use crate::args::CommandArgs;
use crate::errors::AppError;
use crate::fasta::index_fasta;
use crate::utils::{NTS, u8_to_nt};

pub fn bam_parse(args: &CommandArgs) -> Result<(), AppError> {
    if !&args.bam.exists() {
        return Err(AppError::FileDoesNotExistError(
            args.bam.display().to_string(),
        ));
    }

    info!("Reading bam file...");
    let mut bam = bam::Reader::from_path(&args.bam)?;
    bam.set_threads(args.threads)
        .expect("Failed to set num threads for htlib bam reader.");

    let mut pileups = bam.pileup();

    info!("Parsing indexed fasta file...");
    let faidx_reader = faidx::Reader::from_path(&args.fasta)?;
    let seq_names: Vec<String> = faidx_reader.seq_names()?;
    let fasta_index = index_fasta(&faidx_reader, &seq_names);

    // Create writer - either to file or stdout
    let mut writer: BufWriter<Box<dyn Write>> = match &args.output {
        Some(path) => {
            let file = File::create(path).expect("Failed to create output file");
            BufWriter::new(Box::new(file))
        }
        None => BufWriter::new(Box::new(io::stdout())),
    };

    writeln!(
        writer,
        "contig_name\tnt\tpos\tnt_context\tdepth\tref_cov\trev_cov_abs\tdeletion_fraction\tA\tT\tC\tG"
    )
    .expect("Failed to write header");

    // Each nucleotide position.
    while let Some(Ok(pileup)) = pileups.next() {
        let depth = pileup.depth() as usize;
        if depth < args.min_read_depth as usize {
            continue;
        }

        // Extract reference nucleotide at this position.
        // NOTE - It is yet unknown what is the most efficient way to extract a nt from a given
        // reference name and position. The current implementation has a pre-constructed hashmap
        // of structure {"name"; seq}, for which we extract the nt per reference name and position.
        let pos = pileup.pos();
        let reference_id = pileup.tid();
        let reference_name = faidx_reader.seq_name(reference_id as i32)?;
        let fasta_info = fasta_index
            .get(reference_name.as_str())
            .ok_or(AppError::MissingFastaID(reference_name.clone()))?;

        // reference nucleotide.
        let ref_nt = &fasta_info.seq[pos as usize];
        let seq_len = fasta_info.length;

        // For the given position, extract flanking bases to provide context.
        let context_start: usize = max(pos as i32 - args.num_flanking_bases as i32, 0) as usize;
        let context_end: usize =
            min(pos + 1 + (args.num_flanking_bases as u32), seq_len as u32) as usize;

        let context = &fasta_info.seq[context_start..context_end];

        // convert to base.
        let ref_nt_char = u8_to_nt(ref_nt).unwrap();

        // Hashmap for nt counts.
        let random_state = RandomState::new();
        let mut map: HashMap<char, usize> =
            HashMap::with_capacity_and_hasher(NTS.len(), random_state);

        // We can probably re-use this instead of re-creating it several times.
        NTS.iter().for_each(|nt| {
            map.insert(*nt, 0_usize);
        });

        // Alignment of a single read against a single position.
        for alignment in pileup.alignments() {
            if alignment.record().is_unmapped() {
                continue;
            }

            // We just increment num deletions and then skip to next alignment
            if alignment.is_del() || alignment.is_refskip() {
                map.entry('-').and_modify(|c| *c += 1).or_insert(1);
                continue;
            };

            // Position in read.
            let read_pos = alignment
                .qpos()
                .expect("unreachable, we have checked for .id_del() and is_refskip().");

            // Apparently, there are records where seq len is zero.
            // Not sure what this about, need to browse htslib reference.
            if alignment.record().seq_len() == 0 {
                continue;
            }

            let base = alignment.record().seq()[read_pos];
            let read_nt = u8_to_nt(&base).unwrap();

            match map.get_mut(&read_nt) {
                Some(x) => *x += 1,
                _ => continue,
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
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                reference_name,
                ref_nt_char,
                pos,
                std::str::from_utf8(context).unwrap(),
                depth,
                ref_frac_cov,
                ref_abs_frac_cov,
                deletion_fraction,
                map.get(&'A').expect("Missing `A`"),
                map.get(&'T').expect("Missing `T`"),
                map.get(&'C').expect("Missing `C`"),
                map.get(&'G').expect("Missing `G`")
            )
            .expect("Failed to write output");
        }
    }

    writer.flush().expect("Failed to flush output");

    Ok(())
}
