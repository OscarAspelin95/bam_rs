use crate::errors::AppError;
use log::info;
use rust_htslib::bam::{Read, Reader};
use rust_htslib::faidx;
use std::path::Path;

/// Extract reference names from the BAM header, ordered by tid.
pub fn get_reference_names(reader: &Reader) -> Vec<String> {
    let header = reader.header();

    (0..header.target_count())
        .map(|i| String::from_utf8_lossy(header.tid2name(i)).into_owned())
        .collect()
}

/// Load all reference sequences from the FASTA, ordered by BAM header tid.
/// Exits with an error if any BAM contig is missing from the FASTA.
pub fn load_reference_sequences(
    fasta: &Path,
    ref_names: &[String],
) -> Result<Vec<Vec<u8>>, AppError> {
    // Build .fai index if it doesn't exist.
    let fai_path = fasta.with_extension(format!(
        "{}.fai",
        fasta.extension().unwrap_or_default().to_string_lossy()
    ));

    if !fai_path.exists() {
        info!("building FASTA index for `{}`...", fasta.display());
        faidx::build(fasta).map_err(|e| AppError::FaidxBuildError(e.to_string()))?;
    }

    let faidx_reader = faidx::Reader::from_path(fasta)?;

    ref_names
        .iter()
        .enumerate()
        .map(|(tid, name)| {
            let len = faidx_reader.fetch_seq_len(name);

            if len == 0 {
                return Err(AppError::ReferenceMismatchError(format!(
                    "contig `{}` (tid={}) from BAM header not found in FASTA `{}`",
                    name,
                    tid,
                    fasta.display()
                )));
            }

            let seq = faidx_reader.fetch_seq(name, 0, (len - 1) as usize)?;
            Ok(seq)
        })
        .collect()
}
