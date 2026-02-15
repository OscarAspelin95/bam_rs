use thiserror::Error;

#[derive(Debug, Error)]
pub enum AppError {
    #[error("file does not exist: {0}")]
    FileDoesNotExistError(String),

    #[error("Htslib error: {0}")]
    HtslibError(#[from] rust_htslib::errors::Error),

    #[error("Invalid nucleotide: {0}")]
    InvalidNtError(u8),

    #[error("Missing FASTA id: {0}")]
    MissingFastaID(String),
}
