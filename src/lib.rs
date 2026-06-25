//! Rust port of MACS3 `callpeak`: single-end, no-control peak calling.
//!
//! Reads 5' tags from BAM, predicts the fragment length by cross-correlation
//! (or takes `--extsize`), piles up extended tags, derives a dynamic local
//! Poisson background, and emits `narrowPeak` + `summits.bed`.

pub mod model;
pub mod peaks;
pub mod pileup;
pub mod tags;

use std::path::PathBuf;

use rsomics_common::{Result, RsomicsError};

/// callpeak options (single-end default path subset).
pub struct CallPeakOpts {
    pub treatment: Vec<PathBuf>,
    pub control: Vec<PathBuf>,
    pub keep_dup: String,
    pub name: String,
    /// Effective genome size.
    pub gsize: f64,
    /// Skip model building and use `extsize` directly.
    pub nomodel: bool,
    pub extsize: i32,
    /// Output directory for `<name>_peaks.narrowPeak` / `<name>_summits.bed`.
    pub outdir: PathBuf,
}

/// Resolve `--keep-dup` to a per-position cap (0 = keep all).
fn resolve_keep_dup(spec: &str) -> Result<u32> {
    match spec {
        "all" => Ok(0),
        "auto" => Err(RsomicsError::InvalidInput(
            "--keep-dup auto is not supported; pass an integer cap or \"all\"".into(),
        )),
        s => s
            .parse::<u32>()
            .map_err(|_| RsomicsError::InvalidInput(format!("invalid --keep-dup: {s}"))),
    }
}

/// Run no-control single-end callpeak: load + dedup tags, find the fragment
/// length, pile up, call peaks, and write the narrowPeak / summits outputs.
pub fn run_callpeak(opts: &CallPeakOpts) -> Result<()> {
    if !opts.control.is_empty() {
        return Err(RsomicsError::InvalidInput(
            "control (-c) peak calling is not supported; omit -c for no-control callpeak".into(),
        ));
    }
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

    let d = if opts.nomodel {
        opts.extsize
    } else {
        let m = model::build(&treat, opts.gsize, 300, 5.0, 50.0, 20)?;
        eprintln!(
            "predicted fragment length d = {} (alternatives {:?})",
            m.d, m.alternative_d
        );
        m.d
    };
    eprintln!("fragment length d = {d}");

    let treat_track = pileup::treat_pileup_raw(&treat, d);
    let control_track = pileup::control_lambda_raw(&treat, d, opts.gsize, 10000);
    let called = peaks::call_peaks(&treat, d, &treat_track, &control_track);
    eprintln!("{} peaks called", called.len());
    peaks::write_outputs(&called, &treat.names, &opts.name, &opts.outdir)?;
    Ok(())
}
