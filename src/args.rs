use std::path::PathBuf;

use clap::Parser;

use crate::errors::AppError;

fn validate_bam(bam: &str) -> Result<PathBuf, AppError> {
    let path = PathBuf::from(bam);

    if !path.exists() {
        return Err(AppError::FileDoesNotExistError(format!(
            "`{}` does not exist",
            path.display()
        )));
    }

    match path.extension().and_then(|p| p.to_str()) {
        Some("bam") => Ok(path),
        _ => Err(AppError::InvalidFileExtensionError(
            path.display().to_string(),
        )),
    }
}

fn validate_fasta(fasta: &str) -> Result<PathBuf, AppError> {
    let path = PathBuf::from(fasta);

    if !path.exists() {
        return Err(AppError::FileDoesNotExistError(format!(
            "`{}` does not exist",
            path.display()
        )));
    }

    match path.extension().and_then(|p| p.to_str()) {
        Some("fa") | Some("fasta") | Some("fna") => Ok(path),
        _ => Err(AppError::InvalidFileExtensionError(
            path.display().to_string(),
        )),
    }
}

#[derive(Debug, Parser)]
pub struct Args {
    #[arg(short, long, value_parser = validate_bam)]
    pub bam: PathBuf,

    #[arg(short, long, value_parser = validate_fasta)]
    pub fasta: PathBuf,

    #[arg(short, long)]
    pub outfile: PathBuf,

    /// Number of reference bases to include on each side of the pileup position in the `ref` column.
    #[arg(short = 'c', long, default_value_t = 0)]
    pub context: usize,

    /// Only output positions where frac_aln_canonical is less than or equal to this value.
    /// Useful for filtering to heterogeneous positions where reads disagree with the reference.
    #[arg(long)]
    pub min_alternate_frac_canonical: Option<f64>,
}
