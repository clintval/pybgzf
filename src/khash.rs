//! The bucket order of htslib's `khash` integer hash table.
//!
//! htslib stores each reference's bins in a `khash` table and writes them to an index in bucket
//! order, so reproducing that order byte for byte means replaying every insertion, including the
//! table's growth, exactly as `kh_put` and `kh_resize` in htslib's `khash.h` perform them.

const UPPER: f64 = 0.77;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Slot {
    Empty,
    Deleted,
    Live,
}

/// Tracks where `khash` would place each key, without storing any values.
#[derive(Default)]
pub struct KhashOrder {
    keys: Vec<u32>,
    slots: Vec<Slot>,
    size: u32,
    occupied: u32,
    upper_bound: u32,
}

fn upper_bound(buckets: u32) -> u32 {
    (f64::from(buckets) * UPPER + 0.5) as u32
}

fn round_up_to_power_of_two(value: u32) -> u32 {
    let mut x = value.wrapping_sub(1);
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x.wrapping_add(1)
}

impl KhashOrder {
    /// Inserts `key` as `kh_put` does, growing or rehashing the table first when it is full.
    pub fn put(&mut self, key: u32) {
        if self.occupied >= self.upper_bound {
            let buckets = self.keys.len() as u32;
            if buckets > self.size << 1 {
                self.resize(buckets - 1);
            } else {
                self.resize(buckets + 1);
            }
        }
        let buckets = self.keys.len();
        let mask = buckets - 1;
        let mut site = buckets;
        let mut x = buckets;
        let mut i = key as usize & mask;
        if self.slots[i] == Slot::Empty {
            x = i;
        } else {
            let last = i;
            let mut step = 0;
            while self.slots[i] != Slot::Empty
                && (self.slots[i] == Slot::Deleted || self.keys[i] != key)
            {
                if self.slots[i] == Slot::Deleted {
                    site = i;
                }
                step += 1;
                i = (i + step) & mask;
                if i == last {
                    x = site;
                    break;
                }
            }
            if x == buckets {
                x = if self.slots[i] == Slot::Empty && site != buckets {
                    site
                } else {
                    i
                };
            }
        }
        match self.slots[x] {
            Slot::Empty => {
                self.keys[x] = key;
                self.slots[x] = Slot::Live;
                self.size += 1;
                self.occupied += 1;
            }
            Slot::Deleted => {
                self.keys[x] = key;
                self.slots[x] = Slot::Live;
                self.size += 1;
            }
            Slot::Live => {}
        }
    }

    fn resize(&mut self, requested: u32) {
        let buckets = round_up_to_power_of_two(requested).max(4);
        if self.size >= upper_bound(buckets) {
            return;
        }
        let (old, new) = (self.keys.len(), buckets as usize);
        if old < new {
            self.keys.resize(new, 0);
        }
        let mask = new - 1;
        let mut slots = vec![Slot::Empty; new];
        for j in 0..old {
            if self.slots[j] != Slot::Live {
                continue;
            }
            let mut key = self.keys[j];
            self.slots[j] = Slot::Deleted;
            loop {
                let mut i = key as usize & mask;
                let mut step = 0;
                while slots[i] != Slot::Empty {
                    step += 1;
                    i = (i + step) & mask;
                }
                slots[i] = Slot::Live;
                if i < old && self.slots[i] == Slot::Live {
                    std::mem::swap(&mut key, &mut self.keys[i]);
                    self.slots[i] = Slot::Deleted;
                } else {
                    self.keys[i] = key;
                    break;
                }
            }
        }
        self.keys.truncate(new);
        self.slots = slots;
        self.occupied = self.size;
        self.upper_bound = upper_bound(buckets);
    }

    /// Returns the keys in the order `kh_begin` to `kh_end` visits them.
    pub fn keys(&self) -> impl Iterator<Item = u32> + '_ {
        self.keys
            .iter()
            .zip(&self.slots)
            .filter(|(_, slot)| **slot == Slot::Live)
            .map(|(key, _)| *key)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn order(keys: &[u32]) -> Vec<u32> {
        let mut table = KhashOrder::default();
        for &key in keys {
            table.put(key);
        }
        table.keys().collect()
    }

    #[test]
    fn small_tables_place_keys_by_their_low_bits() {
        assert_eq!(order(&[4681, 4682, 37450]), vec![4681, 4682, 37450]);
        assert_eq!(order(&[6, 5]), vec![5, 6]);
    }

    #[test]
    fn repeated_keys_are_stored_once() {
        assert_eq!(order(&[9, 9, 9, 1]), vec![9, 1]);
    }

    #[test]
    fn growth_rehashes_every_key() {
        let keys: Vec<u32> = (0..100).map(|i| 4681 + i * 7).collect();
        let mut expected = keys.clone();
        let placed = order(&keys);
        expected.sort_unstable();
        let mut sorted = placed.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, expected);
        assert_eq!(placed.len(), 100);
    }

    #[test]
    fn rounds_up_like_kroundup32() {
        assert_eq!(round_up_to_power_of_two(0), 0);
        assert_eq!(round_up_to_power_of_two(1), 1);
        assert_eq!(round_up_to_power_of_two(5), 8);
        assert_eq!(round_up_to_power_of_two(8), 8);
    }
}
