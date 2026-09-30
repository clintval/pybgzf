//! Deciding the column layout from content while it streams.

use crate::columns::{Columns, Kind};

const SAM_HEADERS: [&[u8]; 5] = [b"@HD", b"@SQ", b"@RG", b"@PG", b"@CO"];

/// Decides a layout from header lines or, failing that, from the first data line.
#[derive(Default)]
pub struct Sniffer {
    first_bytes: Vec<u8>,
    last_track_line: u64,
    bed_only: bool,
}

fn is_sam_header(line: &[u8]) -> bool {
    SAM_HEADERS
        .iter()
        .any(|tag| line.starts_with(tag) && matches!(line.get(3), None | Some(b'\t')))
}

fn integer(field: &[u8]) -> Option<u64> {
    if field.is_empty() || !field.iter().all(u8::is_ascii_digit) {
        return None;
    }
    std::str::from_utf8(field).ok()?.parse().ok()
}

fn looks_like_gff(fields: &[&[u8]]) -> bool {
    fields.len() == 9
        && matches!((integer(fields[3]), integer(fields[4])), (Some(start), Some(end)) if start <= end)
        && matches!(fields[6], b"+" | b"-" | b"." | b"?")
        && matches!(fields[7], b"." | b"0" | b"1" | b"2")
}

fn looks_like_bed(fields: &[&[u8]]) -> bool {
    fields.len() >= 3
        && matches!((integer(fields[1]), integer(fields[2])), (Some(start), Some(end)) if start <= end)
}

fn looks_like_bed2(fields: &[&[u8]]) -> bool {
    fields.len() == 2 && integer(fields[1]).is_some()
}

impl Sniffer {
    /// Returns a sniffer that only chooses between BED and BED2, for files named as BED.
    pub fn bed() -> Self {
        Self {
            bed_only: true,
            ..Self::default()
        }
    }

    /// Looks at the next line, without its line terminator, and returns the layout once one can
    /// be decided. Errors if the line is data that matches no known layout, or if an earlier line
    /// would not be a header line under the decided layout.
    pub fn push(&mut self, line: &[u8]) -> Result<Option<Columns>, String> {
        let line_number = self.first_bytes.len() as u64 + 1;
        let headers = !self.bed_only;
        let decided = if headers && line.starts_with(b"##fileformat=VCF") {
            Columns::vcf()
        } else if headers && is_sam_header(line) {
            Columns::sam()
        } else {
            let columns = if headers && line.starts_with(b"##gff-version") {
                Columns::gff()
            } else if line.starts_with(b"track ") || line.starts_with(b"browser ") {
                self.last_track_line = line_number;
                self.first_bytes.push(line[0]);
                return Ok(None);
            } else if line.first() == Some(&b'#') {
                self.first_bytes.push(b'#');
                return Ok(None);
            } else {
                let fields: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();
                if self.bed_only && fields.len() == 2 {
                    Columns::bed2()
                } else if self.bed_only {
                    Columns::bed()
                } else if looks_like_gff(&fields) {
                    Columns::gff()
                } else if looks_like_bed(&fields) {
                    Columns::bed()
                } else if looks_like_bed2(&fields) {
                    Columns::bed2()
                } else {
                    return Err(format!(
                        "line {line_number} does not look like BED, GFF, SAM, or VCF: {:?}; pass columns explicitly",
                        String::from_utf8_lossy(line)
                    ));
                }
            };
            Columns {
                skip_lines: self.last_track_line,
                ..columns
            }
        };
        for (i, &first) in self.first_bytes.iter().enumerate() {
            let number = i as u64 + 1;
            if number > decided.skip_lines && first != decided.meta_char {
                return Err(format!(
                    "line {number} is not a header line for the {} layout inferred at line {line_number}; pass columns explicitly",
                    describe(&decided)
                ));
            }
        }
        Ok(Some(decided))
    }
}

fn describe(columns: &Columns) -> &'static str {
    match (columns.kind, columns.start, columns.end) {
        (Kind::Vcf, ..) => "VCF",
        (Kind::Sam, ..) => "SAM",
        (_, 4, _) => "GFF",
        (_, _, None) => "BED2",
        _ => "BED",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sniff(lines: &[&str]) -> Result<Columns, String> {
        let mut sniffer = Sniffer::default();
        for line in lines {
            if let Some(columns) = sniffer.push(line.as_bytes())? {
                return Ok(columns);
            }
        }
        Err("undecided".into())
    }

    #[test]
    fn headers_decide_immediately() {
        assert_eq!(sniff(&["##fileformat=VCFv4.3"]), Ok(Columns::vcf()));
        assert_eq!(sniff(&["@HD\tVN:1.6"]), Ok(Columns::sam()));
        assert_eq!(sniff(&["@CO"]), Ok(Columns::sam()));
        assert_eq!(sniff(&["##gff-version 3"]), Ok(Columns::gff()));
    }

    #[test]
    fn data_lines_decide_between_gff_and_bed() {
        let gff = "chr1\tsrc\tgene\t1\t10\t.\t+\t.\tID=a";
        assert_eq!(sniff(&["# comment", gff]), Ok(Columns::gff()));
        let gtf = "chr1\tsrc\texon\t1\t10\t0.5\t-\t0\tgene_id \"a\";";
        assert_eq!(sniff(&[gtf]), Ok(Columns::gff()));
        assert_eq!(sniff(&["chr1\t0\t10\tname"]), Ok(Columns::bed()));
    }

    #[test]
    fn two_column_lines_are_bed2() {
        assert_eq!(sniff(&["#c", "chr1\t5"]), Ok(Columns::bed2()));
        assert!(sniff(&["chr1\tfive"]).is_err());
    }

    #[test]
    fn bed_files_choose_between_bed_and_bed2() {
        let mut sniffer = Sniffer::bed();
        assert_eq!(sniffer.push(b"##fileformat=VCFv4.2"), Ok(None));
        assert_eq!(sniffer.push(b"chr1\t5"), Ok(Some(Columns::bed2())));
        let mut sniffer = Sniffer::bed();
        assert_eq!(sniffer.push(b"chr1\t5\t6\tx"), Ok(Some(Columns::bed())));
    }

    #[test]
    fn track_and_browser_lines_are_skipped() {
        let lines = [
            "browser position chr1:1-100",
            "track name=x",
            "#chrom\tstart\tend",
            "chr1\t0\t1",
        ];
        assert_eq!(
            sniff(&lines),
            Ok(Columns {
                skip_lines: 2,
                ..Columns::bed()
            })
        );
    }

    #[test]
    fn ambiguous_lines_are_errors() {
        assert!(sniff(&["chr1\t5\t2"]).unwrap_err().contains("line 1"));
        assert!(sniff(&["#x", "hello"]).unwrap_err().contains("line 2"));
        assert!(sniff(&[""]).is_err());
    }

    #[test]
    fn earlier_lines_must_be_headers_under_the_decision() {
        assert!(
            sniff(&["track name=x", "##fileformat=VCFv4.2"])
                .unwrap_err()
                .contains("line 1")
        );
        assert!(
            sniff(&["# note", "@HD\tVN:1.6"])
                .unwrap_err()
                .contains("SAM")
        );
    }
}
