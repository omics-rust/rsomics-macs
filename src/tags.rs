//! 5'-end tag tracks loaded from BAM — the MACS3 `FWTrack` equivalent.
//!
//! Per `.autopilot/state/macs-spec.md` §2: the + strand 5' end is the alignment
//! start; the - strand 5' end is the alignment end (start + reference span). A
//! tag is one 5' position; `filter_dup` caps duplicates at a position.

use std::collections::HashMap;
use std::path::Path;

use rsomics_bamio::raw::{self, RawRecord};
use rsomics_common::{Result, RsomicsError};

const FLAG_REVERSE: u16 = 0x0010;

/// MACS3 `bam_fw_binary_parse` SE drop mask = unmapped | secondary | QC-fail |
/// supplementary (0x4 | 0x100 | 0x200 | 0x800 = 2820).
const SKIP_FLAGS: u16 = 2820;

/// True if MACS3 drops this single-end alignment (Parser.py
/// `bam_fw_binary_parse`). For a paired read it additionally drops the second
/// mate (0x80) / mate-unmapped (0x8) / not-proper-pair (missing 0x2).
fn macs_se_skip(flags: u16) -> bool {
    flags & SKIP_FLAGS != 0 || (flags & 1 != 0 && (flags & 136 != 0 || flags & 2 == 0))
}

/// Per-(reference id, strand) sorted 5'-end positions.
#[derive(Default)]
pub struct Tags {
    pub plus: HashMap<i32, Vec<i32>>,
    pub minus: HashMap<i32, Vec<i32>>,
    /// Tag size = `int(mean l_seq of the first 10 records)`, per MACS3 `tsize()`.
    pub tsize: u32,
    /// Total tags after the most recent `filter_dup`.
    pub total: u64,
}

/// Reference span = sum of CIGAR ops consuming the reference (M/D/N/=/X).
fn reference_span(rec: &RawRecord) -> i32 {
    let mut span: i64 = 0;
    for (op, len) in rec.cigar_ops() {
        if matches!(op, 0 | 2 | 3 | 7 | 8) {
            span += i64::from(len);
        }
    }
    span as i32
}

/// Load 5' tags from `path`, skipping unmapped reads.
pub fn load_bam(path: &Path) -> Result<Tags> {
    let mut reader = rsomics_bamio::open_parallel(path)?;
    reader.read_header().map_err(RsomicsError::Io)?;
    let inner = reader.get_mut();
    let mut rec = RawRecord::default();
    let mut tags = Tags::default();
    let mut tsize_sum = 0u64;
    let mut tsize_n = 0u32;
    loop {
        let n = raw::read_record(inner, &mut rec)?;
        if n == 0 {
            break;
        }
        let flags = rec.flags();
        // MACS3 BAM tsize() = int(mean l_seq of the first 10 records in file
        // order, regardless of flag).
        if tsize_n < 10 {
            tsize_sum += rec.sequence_len() as u64;
            tsize_n += 1;
        }
        if macs_se_skip(flags) {
            continue;
        }
        let tid = rec.reference_sequence_id();
        if flags & FLAG_REVERSE != 0 {
            let p = rec.alignment_start() + reference_span(&rec);
            tags.minus.entry(tid).or_default().push(p);
        } else {
            tags.plus
                .entry(tid)
                .or_default()
                .push(rec.alignment_start());
        }
    }
    if tsize_n > 0 {
        tags.tsize = u32::try_from(tsize_sum / u64::from(tsize_n)).unwrap_or(0);
    }
    Ok(tags)
}

impl Tags {
    /// Tag count across both strands before any dedup.
    #[must_use]
    pub fn raw_count(&self) -> u64 {
        self.plus
            .values()
            .chain(self.minus.values())
            .map(|v| v.len() as u64)
            .sum()
    }

    /// Sort each (reference, strand) array and keep at most `keep_dup` tags at
    /// any identical 5' position (`keep_dup == 0` keeps all). Recomputes `total`.
    pub fn filter_dup(&mut self, keep_dup: u32) {
        let mut total = 0u64;
        for arr in self.plus.values_mut().chain(self.minus.values_mut()) {
            arr.sort_unstable();
            cap_runs(arr, keep_dup);
            total += arr.len() as u64;
        }
        self.total = total;
    }
}

/// Keep at most `cap` copies of each run of equal values in a sorted vec.
fn cap_runs(arr: &mut Vec<i32>, cap: u32) {
    if cap == 0 {
        return;
    }
    let mut w = 0usize;
    let mut i = 0usize;
    while i < arr.len() {
        let v = arr[i];
        let mut kept = 0u32;
        while i < arr.len() && arr[i] == v {
            if kept < cap {
                arr[w] = v;
                w += 1;
                kept += 1;
            }
            i += 1;
        }
    }
    arr.truncate(w);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cap_runs_keep_one() {
        let mut a = vec![1, 1, 1, 2, 3, 3];
        cap_runs(&mut a, 1);
        assert_eq!(a, vec![1, 2, 3]);
    }

    #[test]
    fn cap_runs_keep_two_and_all() {
        let mut a = vec![5, 5, 5, 5];
        cap_runs(&mut a, 2);
        assert_eq!(a, vec![5, 5]);
        let mut b = vec![5, 5, 5, 5];
        cap_runs(&mut b, 0);
        assert_eq!(b, vec![5, 5, 5, 5]);
    }
}
