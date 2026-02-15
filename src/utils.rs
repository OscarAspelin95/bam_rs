use crate::errors::AppError;

pub const NTS: [char; 6] = ['A', 'T', 'C', 'G', 'N', '-'];

/// NOTE - we can probably do this more efficient by having a lookup table
/// statically defined. Also, we might want to distinguish between upper and lowercase nts.
#[inline]
pub fn u8_to_nt(x: &u8) -> Result<char, AppError> {
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
        None => Err(AppError::InvalidNtError(*x)),
    }
}
