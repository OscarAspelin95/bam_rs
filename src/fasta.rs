use rust_htslib::faidx::Reader;
use std::collections::HashMap;

use crate::find_homopolymers;

#[derive(Debug, Clone)]
pub struct Fasta<'a> {
    pub length: usize,
    pub seq: Vec<u8>,
    pub homopolymers: Vec<(&'a str, usize, usize, u8)>,
}

pub fn index_fasta<'a>(reader: &Reader, seq_names: &'a Vec<String>) -> HashMap<&'a str, Fasta<'a>> {
    let fasta_index: HashMap<&'a str, Fasta> = seq_names
        .iter()
        .map(|seq_name| {
            let seq_len = reader.fetch_seq_len(&seq_name);
            let seq = reader
                .fetch_seq(&seq_name, 0, (seq_len + 1) as usize)
                .unwrap();
            let homopolymers = find_homopolymers(&seq, seq_len as usize, &seq_name);

            return (
                seq_name.as_str(),
                Fasta {
                    length: seq_len as usize,
                    seq,
                    homopolymers,
                },
            );
        })
        .collect();

    return fasta_index;
}
