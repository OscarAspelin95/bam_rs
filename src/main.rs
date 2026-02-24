mod args;
mod errors;
mod pileup;
mod reference;
mod seq;

use crate::args::Args;
use crate::errors::AppError;
use crate::pileup::parse;
use clap::Parser;
use log::{self, info};
use simple_logger::SimpleLogger;

fn main() -> Result<(), AppError> {
    SimpleLogger::new()
        .init()
        .expect("failed to initialize logger.");

    let args = Args::parse();

    info!("parsing BAM...");
    parse(
        &args.bam,
        &args.fasta,
        &args.outfile,
        args.context,
        args.min_alternate_frac_canonical,
    )?;

    Ok(())
}
