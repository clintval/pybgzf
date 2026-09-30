//! Finding each line's reference, start, and end, following htslib's `tbx_parse1`.

use memchr::{memchr, memmem};

/// How a line's end is found, and the format code written to the index header.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Generic,
    Sam,
    Vcf,
}

impl Kind {
    pub fn name(self) -> &'static str {
        match self {
            Kind::Generic => "generic",
            Kind::Sam => "sam",
            Kind::Vcf => "vcf",
        }
    }
}

/// Which columns hold each line's reference, start, and end, like `tabix -s -b -e -0 -c -S`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Columns {
    pub refname: usize,
    pub start: usize,
    pub end: Option<usize>,
    pub zero_based: bool,
    pub meta_char: u8,
    pub skip_lines: u64,
    pub kind: Kind,
}

/// A line's reference name and its half-open, zero-based interval.
#[derive(Debug, PartialEq, Eq)]
pub struct Interval<'a> {
    pub name: &'a [u8],
    pub beg: i64,
    pub end: i64,
}

impl Columns {
    pub fn bed() -> Self {
        Self {
            refname: 1,
            start: 2,
            end: Some(3),
            zero_based: true,
            meta_char: b'#',
            skip_lines: 0,
            kind: Kind::Generic,
        }
    }

    pub fn bed2() -> Self {
        Self {
            refname: 1,
            start: 2,
            end: None,
            zero_based: true,
            meta_char: b'#',
            skip_lines: 0,
            kind: Kind::Generic,
        }
    }

    pub fn gff() -> Self {
        Self {
            refname: 1,
            start: 4,
            end: Some(5),
            zero_based: false,
            meta_char: b'#',
            skip_lines: 0,
            kind: Kind::Generic,
        }
    }

    pub fn vcf() -> Self {
        Self {
            refname: 1,
            start: 2,
            end: None,
            zero_based: false,
            meta_char: b'#',
            skip_lines: 0,
            kind: Kind::Vcf,
        }
    }

    pub fn sam() -> Self {
        Self {
            refname: 3,
            start: 4,
            end: None,
            zero_based: false,
            meta_char: b'@',
            skip_lines: 0,
            kind: Kind::Sam,
        }
    }

    /// Checks that the columns describe a layout tabix can record.
    pub fn validate(&self) -> Result<(), String> {
        let limit = i32::MAX as usize;
        if !(1..=limit).contains(&self.refname) || !(1..=limit).contains(&self.start) {
            return Err("column numbers are 1-based and must be positive".into());
        }
        if self.end.is_some_and(|end| !(1..=limit).contains(&end)) {
            return Err("column numbers are 1-based and must be positive".into());
        }
        if self.skip_lines > i32::MAX as u64 {
            return Err("skip_lines is too large".into());
        }
        if self.kind != Kind::Generic && (self.end.is_some() || self.zero_based) {
            return Err(format!(
                "{} columns compute their own end and are 1-based",
                self.kind.name().to_uppercase()
            ));
        }
        if !self.meta_char.is_ascii() {
            return Err("meta_char must be a single ASCII character".into());
        }
        Ok(())
    }

    /// Returns true if the line is a header or comment line that is never indexed.
    pub fn is_meta(&self, line_number: u64, line: &[u8]) -> bool {
        line_number <= self.skip_lines || line.first() == Some(&self.meta_char)
    }

    /// Parses a line, without its line terminator, into an interval.
    pub fn parse<'a>(&self, line: &'a [u8]) -> Result<Interval<'a>, String> {
        match self.kind {
            Kind::Generic => self.parse_generic(line),
            Kind::Sam => parse_sam(line),
            Kind::Vcf => parse_vcf(line),
        }
    }

    fn parse_generic<'a>(&self, line: &'a [u8]) -> Result<Interval<'a>, String> {
        let (sc, bc) = (self.refname, self.start);
        let ec = self.end.unwrap_or(bc);
        let last = sc.max(bc).max(ec);
        let mut name = None;
        let mut beg = -1_i64;
        let mut end = -1_i64;
        for (id, field) in fields(line).enumerate().map(|(i, f)| (i + 1, f)) {
            if id == sc {
                name = Some(field);
            } else if id == bc {
                (beg, end) = start_column(field, bc, ec, self.zero_based, end)?;
            } else if id == ec {
                end = parse_integer(field).ok_or_else(|| not_an_integer(ec, field))?;
            }
            if id >= last {
                break;
            }
        }
        finish(name, beg, end, sc.max(bc))
    }
}

fn fields(line: &[u8]) -> impl Iterator<Item = &[u8]> {
    line.split(|&b| b == b'\t')
}

fn not_an_integer(column: usize, field: &[u8]) -> String {
    format!(
        "column {column} is not an integer: {:?}",
        String::from_utf8_lossy(field)
    )
}

fn start_column(
    field: &[u8],
    bc: usize,
    ec: usize,
    zero_based: bool,
    mut end: i64,
) -> Result<(i64, i64), String> {
    let mut beg = parse_integer(field).ok_or_else(|| not_an_integer(bc, field))?;
    if bc <= ec {
        end = beg;
    }
    if !zero_based {
        beg = beg.saturating_sub(1);
    } else if bc <= ec {
        end = end.saturating_add(1);
    }
    Ok((beg.max(0), end.max(1)))
}

fn finish(name: Option<&[u8]>, beg: i64, end: i64, columns: usize) -> Result<Interval<'_>, String> {
    match name {
        Some(name) if beg >= 0 => {
            if end < 0 {
                return Err(format!("the end position {end} is negative"));
            }
            Ok(Interval { name, beg, end })
        }
        _ => Err(format!("expected at least {columns} tab-separated columns")),
    }
}

/// Parses the leading base-10 integer of `bytes`, saturating, after optional whitespace and
/// sign, returning the value and the number of bytes consumed, or `None` without digits.
/// Unlike htslib, which calls `strtoll` in base 0, a leading `0` or `0x` is not octal or hex.
fn integer_prefix(bytes: &[u8]) -> Option<(i64, usize)> {
    let mut i = bytes.iter().take_while(|b| b.is_ascii_whitespace()).count();
    let negative = match bytes.get(i) {
        Some(b'-') => {
            i += 1;
            true
        }
        Some(b'+') => {
            i += 1;
            false
        }
        _ => false,
    };
    let digits = bytes[i..].iter().take_while(|b| b.is_ascii_digit()).count();
    if digits == 0 {
        return None;
    }
    let mut value: i64 = 0;
    for &b in &bytes[i..i + digits] {
        value = value.saturating_mul(10).saturating_add(i64::from(b - b'0'));
    }
    Some((if negative { -value } else { value }, i + digits))
}

fn parse_integer(bytes: &[u8]) -> Option<i64> {
    integer_prefix(bytes).map(|(value, _)| value)
}

fn integer_or_zero(bytes: &[u8]) -> i64 {
    parse_integer(bytes).unwrap_or(0)
}

fn parse_sam(line: &[u8]) -> Result<Interval<'_>, String> {
    let mut name = None;
    let mut beg = -1_i64;
    let mut end = -1_i64;
    for (id, field) in fields(line).enumerate().map(|(i, f)| (i + 1, f)) {
        match id {
            3 => name = Some(field),
            4 => (beg, end) = start_column(field, 4, 0, false, end)?,
            6 => {
                let mut length: i64 = 0;
                let mut at = 0;
                while at < field.len() {
                    let (count, used) = integer_prefix(&field[at..]).unwrap_or((0, 0));
                    let op = field.get(at + used).map_or(0, u8::to_ascii_uppercase);
                    if matches!(op, b'M' | b'D' | b'N') {
                        length = length.saturating_add(count);
                    }
                    at += used + 1;
                }
                end = beg.saturating_add(if length == 0 { 1 } else { length });
                break;
            }
            _ => {}
        }
    }
    finish(name, beg, end, 4)
}

fn svlen_applies(allele: &[u8]) -> bool {
    allele.len() >= 5
        && allele[0] == b'<'
        && matches!(allele[4], b'>' | b':')
        && matches!(&allele[..4], b"<CNV" | b"<DEL" | b"<DUP" | b"<INV")
        && allele[allele.len() - 1] == b'>'
}

fn info_value<'a>(info: &'a [u8], key: &[u8]) -> Option<&'a [u8]> {
    let at = memmem::find(info, key)?;
    if at == 0 {
        return Some(&info[key.len()..]);
    }
    let mut delimited = Vec::with_capacity(key.len() + 1);
    delimited.push(b';');
    delimited.extend_from_slice(key);
    memmem::find(info, &delimited).map(|at| &info[at + delimited.len()..])
}

#[allow(clippy::too_many_lines)]
fn parse_vcf(line: &[u8]) -> Result<Interval<'_>, String> {
    let mut name = None;
    let mut beg = -1_i64;
    let mut end = -1_i64;
    let mut allele_count = 0_usize;
    let mut svlen_alleles: Vec<bool> = Vec::new();
    let mut uses_svlen = false;
    let mut uses_len = false;
    let mut len_index: Option<usize> = None;
    let (mut reflen, mut svlen, mut fmtlen) = (0_i64, 0_i64, 0_i64);
    for (id, field) in fields(line).enumerate().map(|(i, f)| (i + 1, f)) {
        match id {
            1 => name = Some(field),
            2 => (beg, end) = start_column(field, 2, 0, false, end)?,
            4 => {
                if !field.is_empty() {
                    end = beg.saturating_add(field.len() as i64);
                }
                allele_count += 1;
                svlen_alleles.push(false);
                reflen = field.len() as i64;
            }
            5 => {
                for allele in field.split(|&b| b == b',') {
                    if allele_count >= 65536 {
                        break;
                    }
                    allele_count += 1;
                    let applies = svlen_applies(allele);
                    svlen_alleles.push(applies);
                    if applies {
                        uses_svlen = true;
                    } else if allele == b"<*>" || allele == b"<NON_REF>" {
                        uses_len = true;
                    }
                }
            }
            8 => {
                if let Some(value) = info_value(field, b"END=")
                    && value.first() != Some(&b'.')
                {
                    let info_end = integer_or_zero(value);
                    if info_end > beg {
                        end = info_end;
                    }
                }
                let mut next = info_value(field, b"SVLEN=");
                let mut allele = 1;
                while let Some(value) = next {
                    if allele >= allele_count {
                        break;
                    }
                    let length =
                        if uses_svlen && svlen_alleles.get(allele).copied().unwrap_or(false) {
                            integer_or_zero(value).saturating_abs()
                        } else {
                            1
                        };
                    svlen = svlen.max(length);
                    next = memchr(b',', value).map(|comma| &value[comma + 1..]);
                    allele += 1;
                }
                if !uses_len {
                    break;
                }
            }
            9 if uses_len => {
                len_index = field.split(|&b| b == b':').position(|key| key == b"LEN");
                if len_index.is_none() {
                    break;
                }
            }
            _ if id > 9 && uses_len => {
                if let Some(index) = len_index {
                    let length = field
                        .split(|&b| b == b':')
                        .nth(index)
                        .map_or(0, integer_or_zero);
                    fmtlen = fmtlen.max(length);
                }
            }
            _ => {}
        }
    }
    let longest = reflen.max(svlen).max(fmtlen).saturating_add(beg);
    end = end.max(longest);
    finish(name, beg, end, 2)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn interval(columns: &Columns, line: &str) -> (String, i64, i64) {
        let parsed = columns.parse(line.as_bytes()).unwrap();
        (
            String::from_utf8(parsed.name.to_vec()).unwrap(),
            parsed.beg,
            parsed.end,
        )
    }

    fn named(name: &str, beg: i64, end: i64) -> (String, i64, i64) {
        (name.to_string(), beg, end)
    }

    #[test]
    fn bed_is_zero_based_half_open() {
        assert_eq!(
            interval(&Columns::bed(), "chr1\t10\t20\tname"),
            named("chr1", 10, 20)
        );
    }

    #[test]
    fn bed_without_an_end_column_is_a_point() {
        assert_eq!(interval(&Columns::bed(), "chr1\t10"), named("chr1", 10, 11));
    }

    #[test]
    fn gff_is_one_based_closed() {
        let line = "chr2\tsrc\tgene\t100\t200\t.\t+\t.\tID=a";
        assert_eq!(interval(&Columns::gff(), line), named("chr2", 99, 200));
    }

    #[test]
    fn point_columns_cover_one_base() {
        let one_based = Columns {
            end: None,
            ..Columns::gff()
        };
        assert_eq!(interval(&one_based, "a\tb\tc\t5"), named("a", 4, 5));
        let zero_based = Columns {
            end: None,
            zero_based: true,
            ..Columns::bed()
        };
        assert_eq!(interval(&zero_based, "a\t5"), named("a", 5, 6));
    }

    #[test]
    fn positions_below_one_are_clamped() {
        assert_eq!(interval(&Columns::gff(), "c\ts\tt\t0\t0"), named("c", 0, 0));
    }

    #[test]
    fn missing_columns_and_non_integers_are_errors() {
        assert!(
            Columns::bed()
                .parse(b"chr1")
                .unwrap_err()
                .contains("at least 2")
        );
        assert!(
            Columns::bed()
                .parse(b"chr1\tx\t3")
                .unwrap_err()
                .contains("column 2")
        );
        assert!(
            Columns::bed()
                .parse(b"chr1\t1\tx")
                .unwrap_err()
                .contains("column 3")
        );
    }

    #[test]
    fn sam_end_comes_from_the_cigar() {
        let line = "r\t0\tchr1\t100\t60\t10M5D3I4N2S\t*\t0\t0\tA\tI";
        assert_eq!(interval(&Columns::sam(), line), named("chr1", 99, 99 + 19));
        let unmapped = "r\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI";
        assert_eq!(interval(&Columns::sam(), unmapped), named("*", 0, 1));
    }

    #[test]
    fn vcf_end_comes_from_ref_and_info() {
        let vcf = Columns::vcf();
        assert_eq!(
            interval(&vcf, "1\t100\t.\tACG\tA\t.\t.\t."),
            named("1", 99, 102)
        );
        assert_eq!(
            interval(&vcf, "1\t100\t.\tA\t<DEL>\t.\t.\tEND=500"),
            named("1", 99, 500)
        );
        assert_eq!(
            interval(&vcf, "1\t100\t.\tA\t<DEL>\t.\t.\tSVTYPE=DEL;END=400"),
            named("1", 99, 400)
        );
        assert_eq!(
            interval(&vcf, "1\t100\t.\tA\t<DEL>\t.\t.\tEND=50"),
            named("1", 99, 100)
        );
        assert_eq!(
            interval(&vcf, "1\t100\t.\tA\t<DEL>\t.\t.\tSVLEN=-300"),
            named("1", 99, 399)
        );
        assert_eq!(
            interval(&vcf, "1\t100\t.\tA\tT,<DUP>\t.\t.\tSVLEN=.,50"),
            named("1", 99, 149)
        );
        assert_eq!(
            interval(&vcf, "1\t100\t.\tA\t<INS>\t.\t.\tSVLEN=300"),
            named("1", 99, 100)
        );
    }

    #[test]
    fn gvcf_end_comes_from_format_len() {
        let vcf = Columns::vcf();
        let line = "1\t100\t.\tA\t<*>\t.\t.\t.\tGT:LEN\t0/0:40\t0/0:60";
        assert_eq!(interval(&vcf, line), named("1", 99, 159));
        let without_len = "1\t100\t.\tA\t<*>\t.\t.\t.\tGT\t0/0";
        assert_eq!(interval(&vcf, without_len), named("1", 99, 100));
    }

    #[test]
    fn extreme_positions_saturate() {
        let max = i64::MAX;
        let huge = "99999999999999999999";
        let bed2 = Columns::bed2();
        assert_eq!(interval(&bed2, &format!("c\t{huge}")), named("c", max, max));
        let gff = Columns::gff();
        let line = format!("c\ts\tt\t-{huge}\t{huge}");
        assert_eq!(interval(&gff, &line), named("c", 0, max));
        let sam = format!("r\t0\tc\t{huge}\t60\t{huge}M{huge}D\t*\t0\t0\tA\tI");
        assert_eq!(interval(&Columns::sam(), &sam), named("c", max - 1, max));
        let vcf = Columns::vcf();
        let line = format!("c\t{huge}\t.\tACGT\t<DEL>\t.\t.\tSVLEN=-{huge}");
        assert_eq!(interval(&vcf, &line), named("c", max - 1, max));
        let line = format!("c\t5\t.\tA\t<*>\t.\t.\t.\tLEN\t{huge}");
        assert_eq!(interval(&vcf, &line), named("c", 4, max));
    }

    #[test]
    fn integers_parse_like_strtoll() {
        assert_eq!(integer_prefix(b" 42abc"), Some((42, 3)));
        assert_eq!(integer_prefix(b"-7"), Some((-7, 2)));
        assert_eq!(integer_prefix(b"x"), None);
        assert_eq!(integer_prefix(b""), None);
    }

    #[test]
    fn validation() {
        assert!(Columns::bed().validate().is_ok());
        assert!(
            Columns {
                refname: 0,
                ..Columns::bed()
            }
            .validate()
            .is_err()
        );
        assert!(
            Columns {
                end: Some(3),
                ..Columns::vcf()
            }
            .validate()
            .is_err()
        );
        assert!(
            Columns {
                meta_char: 0xff,
                ..Columns::bed()
            }
            .validate()
            .is_err()
        );
    }
}
