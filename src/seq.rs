use lazy_static::lazy_static;

lazy_static! {
    /// This mapping enables very efficient nucleotide counts later on:
    /// * A/C/G/T maps to index 0,1,2,3.
    /// * N maps to index 4.
    /// * Misc (such as deletions) map to index 5.
    pub static ref ASCII_TO_IDX: [usize; 256] = {
        let mut table = [5usize; 256]; // Default to index 5 (Other/Del)

        table[b'A' as usize] = 0;
        table[b'a' as usize] = 0;
        table[b'C' as usize] = 1;
        table[b'c' as usize] = 1;
        table[b'G' as usize] = 2;
        table[b'g' as usize] = 2;
        table[b'T' as usize] = 3;
        table[b't' as usize] = 3;
        table[b'N' as usize] = 4;
        table[b'n' as usize] = 4;

        table
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(b'A', 0)]
    #[case(b'a', 0)]
    #[case(b'C', 1)]
    #[case(b'c', 1)]
    #[case(b'G', 2)]
    #[case(b'g', 2)]
    #[case(b'T', 3)]
    #[case(b't', 3)]
    #[case(b'N', 4)]
    #[case(b'n', 4)]
    fn ascii_to_idx_nucleotides(#[case] base: u8, #[case] expected: usize) {
        assert_eq!(ASCII_TO_IDX[base as usize], expected);
    }

    #[rstest]
    #[case(b'X')]
    #[case(b'Z')]
    #[case(b'-')]
    #[case(b'.')]
    #[case(0)]
    fn ascii_to_idx_non_nucleotides_map_to_other(#[case] byte: u8) {
        assert_eq!(
            ASCII_TO_IDX[byte as usize], 5,
            "non-nucleotide byte {byte} should map to index 5"
        );
    }
}
