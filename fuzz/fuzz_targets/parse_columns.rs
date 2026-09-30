#![no_main]

use _pybgzf::columns::{Columns, leading_digits};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let digits = data.iter().take_while(|b| b.is_ascii_digit()).count();
    let saturated = data[..digits].iter().fold(0_i64, |value, &b| {
        value.saturating_mul(10).saturating_add(i64::from(b - b'0'))
    });
    assert_eq!(
        leading_digits(data),
        (digits > 0).then_some((saturated, digits))
    );
    let Some((&choice, line)) = data.split_first() else {
        return;
    };
    let generic = Columns {
        refname: usize::from(choice & 3) + 1,
        start: usize::from(choice >> 2 & 3) + 1,
        end: (choice & 0x40 != 0).then(|| usize::from(choice >> 4 & 3) + 1),
        zero_based: choice & 0x80 != 0,
        ..Columns::bed()
    };
    let presets = [
        Columns::bed(),
        Columns::bed2(),
        Columns::gff(),
        Columns::vcf(),
        Columns::sam(),
    ];
    for columns in presets.into_iter().chain([generic]) {
        if let Ok(interval) = columns.parse(line) {
            assert!(interval.beg >= 0 && interval.end >= 0, "{columns:?}");
            assert!(!interval.name.contains(&b'\t'), "{columns:?}");
        }
    }
});
