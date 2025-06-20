use simple_logger::SimpleLogger;

mod bam;
use bam::bam_parse;

mod args;
use args::get_args;

mod utils;

mod fasta;

fn main() {
    SimpleLogger::new().init().unwrap();
    let args = get_args();

    bam_parse(&args);
}
