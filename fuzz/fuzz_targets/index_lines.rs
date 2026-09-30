#![no_main]

use _pybgzf::columns::Columns;
use _pybgzf::writer::{IndexFormat, IndexOptions, Writer};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Some((&choice, lines)) = data.split_first() else {
        return;
    };
    let columns = match choice % 6 {
        0 => Some(Columns::bed()),
        1 => Some(Columns::bed2()),
        2 => Some(Columns::gff()),
        3 => Some(Columns::vcf()),
        4 => Some(Columns::sam()),
        _ => None,
    };
    let format = if choice & 0x80 == 0 {
        IndexFormat::Tabix
    } else {
        IndexFormat::Csi {
            min_shift: u32::from(choice >> 3 & 15) + 1,
            depth: None,
        }
    };
    let path = std::env::temp_dir().join(format!("pybgzf-fuzz-{}.idx", std::process::id()));
    let options = IndexOptions {
        format,
        path,
        columns,
        bed_only: false,
    };
    let mut writer = Writer::new(Vec::new(), 0, 1, Some(options)).expect("valid options");
    let _ = writer.write(lines);
    let _ = writer.finish();
});
