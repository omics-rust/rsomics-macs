//! Rust port of MACS3 `callpeak` (work in progress).
//!
//! The value-exact algorithm is specified in `.autopilot/state/macs-spec.md`.
//! Build sequence: A tag read + dedup [current] → B model d → C pileup +
//! dynamic lambda → D Poisson p / BH q → E peak call + narrowPeak.

pub mod tags;

use std::path::PathBuf;

use rsomics_common::{Result, RsomicsError};

/// callpeak options (single-end default path subset).
pub struct CallPeakOpts {
    pub treatment: Vec<PathBuf>,
    pub control: Vec<PathBuf>,
    pub keep_dup: String,
    pub name: String,
}

/// Resolve `--keep-dup` to a per-position cap (0 = keep all).
fn resolve_keep_dup(spec: &str) -> Result<u32> {
    match spec {
        "all" => Ok(0),
        "auto" => Err(RsomicsError::InvalidInput(
            "--keep-dup auto not yet implemented (binomial cal_max_dup_tags pending)".into(),
        )),
        s => s
            .parse::<u32>()
            .map_err(|_| RsomicsError::InvalidInput(format!("invalid --keep-dup: {s}"))),
    }
}

/// Phase A: load and dedup the treatment track, reporting tag counts to stderr
/// for differential comparison against macs3's reported numbers.
pub fn run_callpeak(opts: &CallPeakOpts) -> Result<()> {
    let cap = resolve_keep_dup(&opts.keep_dup)?;
    let first = opts.treatment.first().ok_or_else(|| {
        RsomicsError::InvalidInput("at least one -t treatment file required".into())
    })?;
    let mut treat = tags::load_bam(first)?;
    let before = treat.raw_count();
    treat.filter_dup(cap);
    eprintln!(
        "treatment tags: {before} total, {} after dup-filter (keep-dup={}), tsize {}",
        treat.total, opts.keep_dup, treat.tsize
    );
    Ok(())
}
