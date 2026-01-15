use clap::Parser;

use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    name = "bam_rs",
    version,
    about = "Find heterogenous alignment positions in BAM files.",
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

    #[arg(short, long, default_value_t = 8)]
    pub threads: usize,

    #[arg(
        short,
        long,
        help = "Output file path. If not provided, writes to stdout."
    )]
    pub output: Option<PathBuf>,

    #[arg(
        long,
        default_value_t = 0.1,
        value_parser = validate_het_frequency,
        help = "Only consider heterogenous position with this fraction or more. E.g., 0.1 means including positions where >=10% of the nucleotides are non reference."
    )]
    pub min_het_frequency: f32,
}

fn validate_het_frequency(s: &str) -> Result<f32, String> {
    let value: f32 = s
        .parse()
        .map_err(|_| format!("'{s}' is not a valid number"))?;
    if !(0.0..=1.0).contains(&value) {
        Err("min_het_frequency must be between 0.0 and 1.0".to_string())
    } else {
        Ok(value)
    }
}

pub fn get_args() -> CommandArgs {
    CommandArgs::parse()
}
