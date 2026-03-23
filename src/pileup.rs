use crate::errors::AppError;
use log::warn;
use rust_htslib::bam::{Read, Reader};
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use crate::reference::{get_reference_names, load_reference_sequences};
use crate::seq::ASCII_TO_IDX;

pub fn parse(
    bam: &Path,
    fasta: &Path,
    outfile: &Path,
    context: usize,
    min_alternate_frac_canonical: Option<f64>,
) -> Result<(), AppError> {
    let mut reader = Reader::from_path(bam)?;
    let mut writer = BufWriter::new(File::create(outfile)?);

    match std::thread::available_parallelism() {
        Ok(n_threads) => reader.set_threads(n_threads.get())?,
        Err(err) => {
            warn!("failed to get available threads. Error: `{err}`")
        }
    }

    let ref_names = get_reference_names(&reader);
    let ref_seqs = load_reference_sequences(fasta, &ref_names)?;

    writer
        .write_all("contig\tpos\tdepth\tref\tA\tC\tG\tT\tN\tdel\tfrac_aln_canonical\tfrac_aln_absolute\tfrac_deletion\tbase_phred\n".as_bytes())
        .expect("failed to write header.");

    reader.pileup().for_each(|pileup| {
        let pileup = match pileup {
            Ok(pileup) => pileup,
            Err(_) => return,
        };

        let tid = pileup.tid() as usize;
        let pos = pileup.pos() as usize;
        let ref_name = &ref_names[tid];
        let ref_seq = &ref_seqs[tid];
        let seq_len = ref_seq.len();

        // Reference base (used for alignment fraction counts).
        let ref_base = ref_seq[pos] as usize;

        // Context window: up to `context` bases on each side, clamped to sequence bounds.
        let ctx_start = pos.saturating_sub(context);
        let ctx_end = (pos + context).min(seq_len - 1);
        let ref_context: String = ref_seq[ctx_start..=ctx_end]
            .iter()
            .map(|&b| b as char)
            .collect();

        let mut counts = [0_usize; 6];
        let mut error_prob_sum: f64 = 0.0;
        let mut phred_count: u32 = 0;

        for alignment in pileup.alignments() {
            let record = alignment.record();

            // Skip entire alignment if not mapped primary.
            if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                continue;
            }

            // Skip deletions (but count occurrence).
            if alignment.is_del() {
                counts[5] += 1;
                continue;
            }

            // Skip reference skips (introns / spliced gaps).
            let read_pos = match alignment.qpos() {
                Some(read_pos) => read_pos,
                None => continue,
            };

            // Phred and base for this alignment at the given position.
            let phred = record.qual()[read_pos];
            let base = record.seq()[read_pos] as usize;

            counts[ASCII_TO_IDX[base]] += 1;
            error_prob_sum += 10_f64.powf(-(phred as f64) / 10.0);
            phred_count += 1;
        }

        let mean_phred = if phred_count > 0 {
            let mean_error = error_prob_sum / phred_count as f64;
            (-10.0 * mean_error.log10()).round() as u8
        } else {
            0
        };

        let [count_a, count_c, count_g, count_t, count_n, count_del] = counts;

        //
        let ref_base_count = counts[ASCII_TO_IDX[ref_base]];
        let canonical_count = count_a + count_c + count_g + count_t;

        // fraction read bases that are the same as the reference base.
        let mut frac_aln_canonical = 0.0;
        if canonical_count > 0 {
            frac_aln_canonical = ref_base_count as f64 / canonical_count as f64;
        }

        // fraction of depth that is the same as the reference base.
        let mut frac_aln_absolute = 0.0;
        if pileup.depth() > 0 {
            frac_aln_absolute = ref_base_count as f64 / pileup.depth() as f64
        }

        // total deletion rate
        let mut frac_deletion = 0.0;
        if pileup.depth() > 0 {
            frac_deletion = count_del as f64 / pileup.depth() as f64;
        }

        if let Some(threshold) = min_alternate_frac_canonical
            && threshold < frac_aln_canonical
        {
            return;
        }

        writer
            .write_all(
                format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{}\n",
                    ref_name,
                    pos,
                    pileup.depth(),
                    ref_context,
                    count_a,
                    count_c,
                    count_g,
                    count_t,
                    count_n,
                    count_del,
                    frac_aln_canonical,
                    frac_aln_absolute,
                    frac_deletion,
                    mean_phred
                )
                .as_bytes(),
            )
            .expect("failed to write results.");
    });

    writer.flush()?;
    Ok(())
}
