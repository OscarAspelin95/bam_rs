use simple_logger::SimpleLogger;

mod bam;
use bam::bam_parse;

mod args;
use args::get_args;

use crate::errors::AppError;

mod errors;
mod fasta;
mod utils;

fn main() -> Result<(), AppError> {
    SimpleLogger::new().init().unwrap();
    let args = get_args();

    bam_parse(&args)?;

    Ok(())
}
