use thiserror::Error;

#[derive(Debug, Error)]
pub enum AppError {
    #[error("File does not exist: `{0}`")]
    FileDoesNotExistError(String),

    #[error("Invalid file extension: `{0}`")]
    InvalidFileExtensionError(String),

    #[error("Reference mismatch: {0}")]
    ReferenceMismatchError(String),

    #[error("HtsLib error: {0}")]
    HtsLibError(#[from] rust_htslib::errors::Error),

    #[error("FASTA index build error: {0}")]
    FaidxBuildError(String),

    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
}
