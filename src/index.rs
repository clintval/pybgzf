//! A binning index built record by record, following htslib's `hts_idx_push` and
//! `hts_idx_finish`, then written with noodles as a tabix or CSI index.

use std::collections::HashMap;
use std::io::{self, Write};

use bstr::BString;
use indexmap::{IndexMap, IndexSet};
use noodles_csi::binning_index::index::header::format::CoordinateSystem;
use noodles_csi::binning_index::index::header::{Format, Header};
use noodles_csi::binning_index::index::reference_sequence::Bin;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles_csi::binning_index::index::reference_sequence::index::{
    BinnedIndex, Index as ReferenceIndex, LinearIndex,
};
use noodles_csi::binning_index::index::{Index, ReferenceSequence};

use crate::columns::{Columns, Kind};
use crate::khash::KhashOrder;

const MIN_MARKER_DISTANCE: u64 = 0x10000;
const UNSET: u64 = u64::MAX;

/// The largest position a binning index with these parameters can hold, saturating at
/// `i64::MAX`.
pub fn max_position(min_shift: u32, depth: u32) -> i64 {
    match min_shift.saturating_add(depth.saturating_mul(3)) {
        shift @ 0..63 => 1_i64 << shift,
        _ => i64::MAX,
    }
}

fn bin_first(level: u32) -> u32 {
    ((1_u32 << (3 * level)) - 1) / 7
}

fn bin_parent(bin: u32) -> u32 {
    (bin - 1) >> 3
}

fn bin_level(bin: u32) -> u32 {
    let (mut level, mut b) = (0, bin);
    while b > 0 {
        level += 1;
        b = bin_parent(b);
    }
    level
}

fn bin_bottom(bin: u32, depth: u32) -> usize {
    let level = bin_level(bin);
    ((bin - bin_first(level)) as usize) << ((depth - level) * 3)
}

/// The bin that holds `[beg, end)`, as in htslib's `hts_reg2bin`.
pub fn reg2bin(beg: i64, end: i64, min_shift: u32, depth: u32) -> u32 {
    let end = end - 1;
    let mut shift = min_shift;
    let mut offset = bin_first(depth);
    for level in (1..=depth).rev() {
        if beg >> shift == end >> shift {
            return offset + (beg >> shift) as u32;
        }
        shift += 3;
        offset -= 1 << (3 * (level - 1));
    }
    0
}

#[derive(Default)]
struct Reference {
    bins: HashMap<u32, Vec<(u64, u64)>>,
    order: KhashOrder,
    linear: Vec<u64>,
    metadata: Option<(u64, u64, u64)>,
}

impl Reference {
    fn add_chunk(&mut self, bin: u32, start: u64, end: u64) {
        self.order.put(bin);
        self.bins.entry(bin).or_default().push((start, end));
    }

    fn set_metadata(&mut self, depth: u32, start: u64, end: u64, count: u64) {
        let id = metadata_bin(depth);
        // htslib puts this key twice, and khash may resize on a put even when the key exists.
        self.order.put(id);
        self.order.put(id);
        self.metadata = Some((start, end, count));
    }

    fn add_to_linear_index(&mut self, beg: i64, end: i64, offset: u64, min_shift: u32) {
        let first = (beg >> min_shift) as usize;
        let last = ((end - 1) >> min_shift) as usize;
        if self.linear.len() < last + 1 {
            self.linear.resize(last + 1, UNSET);
        }
        for slot in &mut self.linear[first..=last] {
            if *slot == UNSET {
                *slot = offset;
            }
        }
    }

    fn backfill_linear_index(&mut self) {
        for i in (0..self.linear.len().saturating_sub(1)).rev() {
            if self.linear[i] == UNSET {
                self.linear[i] = self.linear[i + 1];
            }
        }
    }

    fn compress_bins(&mut self, depth: u32) {
        let bin_count = bin_first(depth + 1);
        for level in (1..=depth).rev() {
            let first = bin_first(level);
            let mut ids: Vec<u32> = self
                .bins
                .keys()
                .copied()
                .filter(|&id| id >= first && id < bin_count)
                .collect();
            ids.sort_unstable();
            for id in ids {
                let Some(chunks) = self.bins.get_mut(&id) else {
                    continue;
                };
                if level < depth {
                    chunks.sort_unstable();
                }
                let span = (chunks[chunks.len() - 1].1 >> 16) - (chunks[0].0 >> 16);
                if span >= MIN_MARKER_DISTANCE || !self.bins.contains_key(&bin_parent(id)) {
                    continue;
                }
                let moved = self.bins.remove(&id).unwrap_or_default();
                self.bins.entry(bin_parent(id)).or_default().extend(moved);
            }
        }
        if let Some(chunks) = self.bins.get_mut(&0) {
            chunks.sort_unstable();
        }
        for chunks in self.bins.values_mut() {
            let mut merged: Vec<(u64, u64)> = Vec::with_capacity(chunks.len());
            for &(start, end) in chunks.iter() {
                match merged.last_mut() {
                    Some(last) if last.1 >> 16 >= start >> 16 => last.1 = last.1.max(end),
                    _ => merged.push((start, end)),
                }
            }
            *chunks = merged;
        }
    }
}

/// Accumulates records, in sorted order, into a binning index.
pub struct IndexBuilder {
    min_shift: u32,
    depth: u32,
    references: Vec<Reference>,
    save_tid: Option<usize>,
    last_tid: Option<usize>,
    save_bin: Option<u32>,
    last_bin: Option<u32>,
    save_offset: u64,
    last_offset: u64,
    first_offset: u64,
    record_count: u64,
}

impl IndexBuilder {
    /// Starts an index whose first record begins at virtual position `first_offset`.
    pub fn new(min_shift: u32, depth: u32, first_offset: u64) -> Self {
        Self {
            min_shift,
            depth,
            references: Vec::new(),
            save_tid: None,
            last_tid: None,
            save_bin: None,
            last_bin: None,
            save_offset: first_offset,
            last_offset: first_offset,
            first_offset,
            record_count: 0,
        }
    }

    /// Adds a record on reference `tid` covering `[beg, end)` that ends at virtual position
    /// `end_offset`. Records must arrive sorted and each reference must be contiguous.
    pub fn push(&mut self, tid: usize, beg: i64, end: i64, end_offset: u64) {
        if self.references.len() < tid + 1 {
            self.references.resize_with(tid + 1, Reference::default);
        }
        if self.last_tid != Some(tid) {
            self.last_tid = Some(tid);
            self.last_bin = None;
        }
        let beg = beg.max(0);
        let end = if end <= 0 { 1 } else { end };
        self.references[tid].add_to_linear_index(beg, end, self.last_offset, self.min_shift);
        let bin = reg2bin(beg, end, self.min_shift, self.depth);
        if self.last_bin != Some(bin) {
            if let (Some(save_tid), Some(save_bin)) = (self.save_tid, self.save_bin) {
                self.references[save_tid].add_chunk(save_bin, self.save_offset, self.last_offset);
                if self.last_bin.is_none() {
                    self.references[save_tid].set_metadata(
                        self.depth,
                        self.first_offset,
                        self.last_offset,
                        self.record_count,
                    );
                    self.record_count = 0;
                    self.first_offset = self.last_offset;
                }
            }
            self.save_offset = self.last_offset;
            self.save_bin = Some(bin);
            self.last_bin = Some(bin);
            self.save_tid = Some(tid);
        }
        self.record_count += 1;
        self.last_offset = end_offset;
    }

    /// Closes the last chunk at `final_offset` and returns the per-reference bins, linear
    /// indexes, and metadata.
    fn finish(mut self, final_offset: u64) -> (u32, u32, Vec<Reference>) {
        if let (Some(tid), Some(bin)) = (self.save_tid, self.save_bin) {
            self.references[tid].add_chunk(bin, self.save_offset, final_offset);
            self.references[tid].set_metadata(
                self.depth,
                self.first_offset,
                final_offset,
                self.record_count,
            );
        }
        for reference in &mut self.references {
            reference.backfill_linear_index();
            reference.compress_bins(self.depth);
        }
        (self.min_shift, self.depth, self.references)
    }
}

fn header(columns: &Columns, names: &IndexSet<Vec<u8>>) -> Header {
    let format = match columns.kind {
        Kind::Sam => Format::Sam,
        Kind::Vcf => Format::Vcf,
        Kind::Generic if columns.zero_based => Format::Generic(CoordinateSystem::Bed),
        Kind::Generic => Format::Generic(CoordinateSystem::Gff),
    };
    noodles_csi::binning_index::index::header::Builder::default()
        .set_format(format)
        .set_reference_sequence_name_index(columns.refname - 1)
        .set_start_position_index(columns.start - 1)
        .set_end_position_index(columns.end.map(|end| end - 1))
        .set_line_comment_prefix(columns.meta_char)
        .set_line_skip_count(columns.skip_lines as u32)
        .set_reference_sequence_names(
            names
                .iter()
                .map(|name| BString::from(name.clone()))
                .collect(),
        )
        .build()
}

fn chunks(pairs: &[(u64, u64)]) -> Vec<Chunk> {
    pairs
        .iter()
        .map(|&(start, end)| Chunk::new(start.into(), end.into()))
        .collect()
}

fn metadata_bin(depth: u32) -> u32 {
    bin_first(depth + 1) + 1
}

/// Returns the bins, and the metadata pseudo-bin, in the order htslib writes them.
fn ordered_bins(reference: &Reference, depth: u32) -> impl Iterator<Item = (u32, Bin)> + '_ {
    let metadata_id = metadata_bin(depth);
    reference.order.keys().filter_map(move |id| {
        let pairs = match (id == metadata_id, reference.metadata) {
            (true, Some((start, end, count))) => &[(start, end), (count, 0)][..],
            _ => reference.bins.get(&id)?,
        };
        Some((id, Bin::new(chunks(pairs))))
    })
}

fn build<I: ReferenceIndex>(
    (min_shift, depth): (u32, u32),
    header: Header,
    references: Vec<ReferenceSequence<I>>,
) -> Index<I> {
    Index::builder()
        .set_min_shift(min_shift as u8)
        .set_depth(depth as u8)
        .set_header(header)
        .set_reference_sequences(references)
        .set_unplaced_unmapped_record_count(0)
        .build()
}

/// Writes the index of `builder`'s records, as CSI when `csi` and tabix otherwise, as BGZF to
/// `sink`.
pub fn write<W: Write>(
    sink: W,
    builder: IndexBuilder,
    final_offset: u64,
    columns: &Columns,
    names: &IndexSet<Vec<u8>>,
    csi: bool,
) -> io::Result<W> {
    let (min_shift, depth, references) = builder.finish(final_offset);
    let header = header(columns, names);
    if !csi {
        let references = references
            .iter()
            .map(|reference| {
                let bins: IndexMap<usize, Bin> = ordered_bins(reference, depth)
                    .map(|(id, bin)| (id as usize, bin))
                    .collect();
                let linear: LinearIndex = reference.linear.iter().map(|&o| o.into()).collect();
                ReferenceSequence::new(bins, linear, None)
            })
            .collect();
        let mut writer = noodles_tabix::io::Writer::new(sink);
        writer.write_index(&build((min_shift, depth), header, references))?;
        return writer.into_inner().finish();
    }
    let bin_count = bin_first(depth + 1);
    let references = references
        .iter()
        .map(|reference| {
            let mut bins = IndexMap::new();
            let mut offsets = BinnedIndex::new();
            for (id, bin) in ordered_bins(reference, depth) {
                let offset = if id < bin_count {
                    let bottom = bin_bottom(id, depth);
                    reference.linear.get(bottom).copied().unwrap_or(0)
                } else {
                    0
                };
                bins.insert(id as usize, bin);
                offsets.insert(id as usize, offset.into());
            }
            ReferenceSequence::new(bins, offsets, None)
        })
        .collect();
    let mut writer = noodles_csi::io::Writer::new(sink);
    writer.write_index(&build((min_shift, depth), header, references))?;
    writer.into_inner().finish()
}

#[cfg(test)]
fn sorted_bins(reference: &Reference) -> Vec<(u32, &Vec<(u64, u64)>)> {
    let mut bins: Vec<_> = reference
        .bins
        .iter()
        .map(|(&id, pairs)| (id, pairs))
        .collect();
    bins.sort_unstable_by_key(|&(id, _)| id);
    bins
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reg2bin_matches_the_specification() {
        assert_eq!(reg2bin(0, 1, 14, 5), 4681);
        assert_eq!(reg2bin(0, 1 << 14, 14, 5), 4681);
        assert_eq!(reg2bin(0, (1 << 14) + 1, 14, 5), 585);
        assert_eq!(reg2bin(0, 1 << 29, 14, 5), 0);
        assert_eq!(reg2bin(1 << 26, (1 << 26) + 1, 14, 5), 4681 + (1 << 12));
    }

    #[test]
    fn bottom_bins() {
        assert_eq!(bin_bottom(0, 5), 0);
        assert_eq!(bin_bottom(1, 5), 0);
        assert_eq!(bin_bottom(2, 5), 1 << 12);
        assert_eq!(bin_bottom(4681, 5), 0);
        assert_eq!(bin_bottom(4682, 5), 1);
    }

    #[test]
    fn records_in_one_bin_share_a_chunk() {
        let mut builder = IndexBuilder::new(14, 5, 0);
        builder.push(0, 0, 10, 20);
        builder.push(0, 5, 15, 40);
        let (_, _, references) = builder.finish(60);
        assert_eq!(references[0].bins[&4681], vec![(0, 60)]);
        assert_eq!(references[0].linear, vec![0]);
        assert_eq!(references[0].metadata, Some((0, 60, 2)));
    }

    #[test]
    fn metadata_is_split_by_reference() {
        let mut builder = IndexBuilder::new(14, 5, 7);
        builder.push(0, 0, 10, 20);
        builder.push(1, 0, 10, 40);
        let (_, _, references) = builder.finish(60);
        assert_eq!(references[0].metadata, Some((7, 20, 1)));
        assert_eq!(references[1].metadata, Some((20, 60, 1)));
    }

    #[test]
    fn linear_index_is_backfilled() {
        let mut builder = IndexBuilder::new(14, 5, 0);
        builder.push(0, 3 << 14, (3 << 14) + 1, 100);
        let (_, _, references) = builder.finish(200);
        assert_eq!(references[0].linear, vec![0, 0, 0, 0]);
    }

    #[test]
    fn small_bins_merge_into_existing_parents() {
        let mut builder = IndexBuilder::new(14, 5, 0);
        builder.push(0, 0, 1 << 15, 10);
        builder.push(0, 1, 2, 20);
        let (_, _, references) = builder.finish(30);
        let ids: Vec<_> = sorted_bins(&references[0])
            .into_iter()
            .map(|(id, _)| id)
            .collect();
        assert_eq!(ids, vec![585]);
        assert_eq!(references[0].bins[&585], vec![(0, 30)]);
    }
}
