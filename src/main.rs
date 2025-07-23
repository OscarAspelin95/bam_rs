use simple_logger::SimpleLogger;

mod bam;
use bam::bam_parse;

mod args;
use args::get_args;

mod fasta;
mod utils;

/// TODO - change to output file instead of writing to stdout.
fn main() {
    SimpleLogger::new().init().unwrap();
    let args = get_args();

    bam_parse(&args);
}
