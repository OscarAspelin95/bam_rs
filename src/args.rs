use clap::Parser;

use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    long_about = "Parses a BAM file of aligned Nanopore reads, outputting heterogenous positions."
)]
pub struct CommandArgs {
    #[arg(short, long, help = "Path to BAM file.")]
    pub bam: PathBuf,

    #[arg(short, long, help = "Path to FASTA file.")]
    pub fasta: PathBuf,

    #[arg(long, default_value = "2", value_parser= clap::value_parser!(u8).range(1..10), help="Extract this many reference nucleotides to the left and right of the heterogenous position")]
    pub num_flanking_bases: u8,

    #[arg(long, default_value = "10",  value_parser= clap::value_parser!(u8).range(1..), help = "Only consider positions with at least this depth.")]
    pub min_read_depth: u8,

    #[arg(
        long,
        default_value_t = 0.1,
        help = "Only consider heterogenous position with this fraction or more. E.g., 0.1 means including positions where >=10% of the nucleotides are non reference."
    )]
    pub min_het_frequency: f32,
}

pub fn get_args() -> CommandArgs {
    let args = CommandArgs::parse();

    if args.min_het_frequency < 0.0 || args.min_het_frequency > 1.0 {
        panic!("Min het frequency must be between 0.0 and 1.0");
    }

    return args;
}
