# bam_rs
Find heterogenous alignment positions in a BAM file (suitable for Nanopore alignments).

## Requirements
- Linux OS (Ubuntu 24.04.2)
- Rust >= 1.88.0

## Install
The easiest way to get started is to download a precompiled Linux binary from the latest [release](https://github.com/OscarAspelin95/bam_rs/releases).

## Install from source
Clone the repository or download the source code. Enter the bam_rs directory and run:

```bash
cargo build --release
```

The generated binary is available in `target/release/bam_rs`.

## Usage
Run with:

```bash
bam_rs --bam <aln.bam> --fasta <sequences.fasta>
```

### Optional arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--output` | stdout | Output file path. If not provided, writes to stdout. |
| `--num-flanking-bases` | 2 | How many reference bases to extract around the heterogenous position. |
| `--min-read-depth` | 10 | Min depth for candidate position. |
| `--threads` | 8 | Num threads to use for Rust HTSlib. |
| `--min-het-frequency` | 0.1 | Minimum heterogenous frequency for candidate position. E.g., 0.1 means consider positions where >= 10% of aligned bases in a position are heterogenous. |

## TODO
- [ ] Add tests.
- [ ] Refactor.
