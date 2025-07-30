# bam_rs
Script for finding heterogenous alignment positions in a BAM file (suitable for Nanopore alignments).

## Requirements
- Linux OS (Ubuntu 24.04.2)
- Rust >= 1.88.0

## Installation
Clone the repository or download the source code. Enter the bam_rs directory and run:<br>
`cargo build --release`

The generated binary is available in `target/release/bam_rs`.

## Usage
Run with:<br>
`bam_rs --bam <aln.bam> --fasta <sequences.fasta>`
