# bam_rs
Find heterogenous alignment positions in a BAM file (suitable for Nanopore alignments).

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

Optional arguments:
<pre>
<b>--num-flanking-bases</b> [2] - How many reference bases to extract around the heterogenous position.

<b>--min_read_depth</b> [10] - Min depth for candidate position.

<b>--threads</b> [8] - Num threads to use for Rust HTSlib.

<b>--min-het-frequency</b> [0.1] - Minimum heterogenous frequency for candidate position. E.g., 0.1 means consider positions where >= 10% of aligned bases in a position are heterogenous.
</pre>
