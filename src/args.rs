

use clap::Parser;

use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(long_about = "Parses a BAM file of aligned Nanopore reads, outputting heterogenous positions.")]
pub struct CommandArgs {
    #[arg(short, long)]
    pub bam: PathBuf,

    #[arg(short, long)]
    pub fasta: PathBuf,

    #[arg(long, default_value = "2", value_parser= clap::value_parser!(u8).range(1..10))]
    pub num_flanking_bases: u8,

    #[arg(long, default_value = "5",  value_parser= clap::value_parser!(u8).range(5..20))]
    pub min_hp_length: u8,

    #[arg(long, default_value = "10",  value_parser= clap::value_parser!(u8).range(1..))]
    pub min_read_depth: u8,

    #[arg(long, default_value_t = 0.9)]
    pub min_het_frequency: f32,

    #[arg(short, long, required = true)]
    pub outfile: PathBuf,
}


pub fn get_args() -> CommandArgs{
    let args = CommandArgs::parse();

    // For now, have extra check for min het frequency.
    if args.min_het_frequency > 1.0{
        panic!("Min het frequency cannot be larger than 1.0");
    }

    return args
}