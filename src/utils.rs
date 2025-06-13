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
