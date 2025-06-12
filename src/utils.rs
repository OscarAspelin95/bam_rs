use thiserror::Error;

#[derive(Error, Debug)]
pub enum BamParseError {
    #[error("Invalid nucleotide")]
    InvalidNtError(u8),
}

/// NOTE - we can probably do this more efficient by having a lookup table
/// statically defined. Also, we might want to distinguish between upper and lowercase nts.
#[inline]
pub fn u8_to_nt(x: &u8) -> Result<char, BamParseError> {
    let nt: Option<char> = match x {
        // A/a.
        b'A' | b'a' => Some('A'),
        // T/t.
        b'T' | b't' => Some('T'),
        // C/c.
        b'C' | b'c' => Some('C'),
        // G/g.
        b'G' | b'g' => Some('G'),
        // N only.
        b'N' => Some('N'),
        _ => None,
    };

    match nt {
        Some(nucleotide) => Ok(nucleotide),
        None => Err(BamParseError::InvalidNtError(*x)),
    }
}

/// Function returns a vector of tuples, where each tuple is on the form (start, end, hp_nt)
/// Where:
/// * hp_nt is in binary format.
/// * start is zero indexed.
/// * end is exclusive.
/// NOTE - not used currently.
pub fn find_homopolymers<'a>(
    nt_str: &[u8],
    length: usize,
    ref_name: &'a str,
) -> Vec<(&'a str, usize, usize, u8)> {
    let mut hp_vec: Vec<(&'a str, usize, usize, u8)> = vec![];

    if nt_str.len() < length {
        return hp_vec;
    }

    let mut i = 0;
    let mut j = 1;

    while i <= nt_str.len() {
        while j < nt_str.len() && nt_str[j] == nt_str[i] {
            j += 1;
        }

        if j - i >= length {
            hp_vec.push((ref_name, i, j, nt_str[i]));
        }

        i = j;
        j += 1;
    }

    return hp_vec;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_homopolymers() {
        // Basic case with a single homopolymer.
        assert_eq!(
            find_homopolymers(b"AAAA", 4, "some_contig"),
            vec![("some_contig", 0, 4, b'A')]
        );

        // Make sure we can catch hps of different nts.
        // We skip T since the length is too short.
        assert_eq!(
            find_homopolymers(b"AGTTTTACGGGGGGCACCCCCCCCCCCCCCCC", 5, "some_contig"),
            vec![("some_contig", 8, 14, b'G'), ("some_contig", 16, 32, b'C')]
        );

        // Make sure we skip short sequences.
        assert_eq!(find_homopolymers(b"ATCG", 3, "some_contig"), vec![]);
    }
}
