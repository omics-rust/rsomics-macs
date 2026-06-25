//! Poisson p-score, BH q-value table, and peak calling — the final MACS3
//! `callpeak` stage.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

use crate::pileup::RawTrack;
use crate::tags::Tags;

/// A combined segment: `(end position, treatment value, control value)`.
type Seg = (i32, f32, f32);

/// A called peak (0-based, half-open `[start, end)`).
pub struct Peak {
    pub tid: i32,
    pub start: i32,
    pub end: i32,
    pub summit: i32,
    pub pileup: f32,
    pub pscore: f32,
    pub qscore: f32,
    pub fc: f32,
}

/// Round half to even at `places` decimals (Python `round`).
fn round_half_even(x: f64, places: i32) -> f64 {
    let factor = 10f64.powi(places);
    let shifted = x * factor;
    let floor = shifted.floor();
    let frac = shifted - floor;
    let rounded = if (frac - 0.5).abs() < 1e-9 {
        if (floor as i64) % 2 == 0 {
            floor
        } else {
            floor + 1.0
        }
    } else {
        shifted.round()
    };
    rounded / factor
}

fn logspace_add(a: f64, b: f64) -> f64 {
    if a > b {
        a + (b - a).exp().ln_1p()
    } else {
        b + (a - b).exp().ln_1p()
    }
}

/// `log10(P(X > k))` for `X ~ Poisson(lbd)`, returned NEGATED (= -log10 tail),
/// rounded to 5 decimals (banker's), matching MACS3 `log10_poisson_cdf_Q`.
fn neg_log10_poisson_q(k: u32, lbd: f64) -> f64 {
    let ln_lbd = lbd.ln();
    let m0 = k + 1;
    let sum_ln_m: f64 = (1..=u64::from(m0)).map(|i| (i as f64).ln()).sum();
    let mut logx = f64::from(m0) * ln_lbd - sum_ln_m;
    let mut residue = logx;
    let mut m = u64::from(m0);
    loop {
        m += 1;
        let logy = logx + ln_lbd - (m as f64).ln();
        let pre = residue;
        residue = logspace_add(pre, logy);
        if (pre - residue).abs() < 1e-5 {
            break;
        }
        logx = logy;
    }
    -round_half_even((residue - lbd) / std::f64::consts::LN_10, 5)
}

/// p-score = `-log10(P(X >= int(treat) | Poisson(ctrl)))`, narrowed to f32 and
/// cached by `(int(treat), ctrl_bits)` as MACS3 does.
fn pscore(treat: f32, ctrl: f32, cache: &mut HashMap<(i32, u32), f32>) -> f32 {
    let ti = treat as i32;
    let key = (ti, ctrl.to_bits());
    if let Some(&v) = cache.get(&key) {
        return v;
    }
    let v = neg_log10_poisson_q(ti.max(0) as u32, f64::from(ctrl)) as f32;
    cache.insert(key, v);
    v
}

/// Combine the raw treatment + control `(end_pos, value)` arrays into
/// `(end, treat, ctrl)` segments (MACS3 `over_two_pv_array`: emit the smaller
/// endpoint, no trailing drain). Uses RAW granularity (not the merged bedGraph)
/// because the summit is the midpoint of the max-treatment SEGMENT.
fn combine(tp: &[i32], tv: &[f32], cp: &[i32], cv: &[f32]) -> Vec<Seg> {
    let mut out = Vec::new();
    let (mut i, mut j) = (0usize, 0usize);
    while i < tp.len() && j < cp.len() {
        let (te, ce) = (tp[i], cp[j]);
        out.push((te.min(ce), tv[i], cv[j]));
        if te < ce {
            i += 1;
        } else if ce < te {
            j += 1;
        } else {
            i += 1;
            j += 1;
        }
    }
    out
}

/// BH q-value table: `pscore(f32 bits) -> qscore`, ranking by base-pair length.
fn build_pqtable(combined: &[&[Seg]], cache: &mut HashMap<(i32, u32), f32>) -> HashMap<u32, f32> {
    let mut bp_by_pscore: HashMap<u32, i64> = HashMap::new();
    for segs in combined {
        let mut pre = 0i32;
        for &(end, tv, cv) in *segs {
            let ps = pscore(tv, cv, cache);
            *bp_by_pscore.entry(ps.to_bits()).or_insert(0) += i64::from(end - pre);
            pre = end;
        }
    }
    let n: i64 = bp_by_pscore.values().sum();
    let mut uniq: Vec<f32> = bp_by_pscore.keys().map(|&b| f32::from_bits(b)).collect();
    uniq.sort_unstable_by(|a, b| b.partial_cmp(a).unwrap()); // descending
    let f = -(n as f64).log10();
    let mut pqtable: HashMap<u32, f32> = HashMap::new();
    let mut k: i64 = 1;
    let mut pre_q = f64::from(i32::MAX);
    let mut broke = uniq.len();
    for (idx, &v) in uniq.iter().enumerate() {
        let l = bp_by_pscore[&v.to_bits()];
        let mut q = f64::from(v) + ((k as f64).log10() + f);
        if q > pre_q {
            q = pre_q;
        }
        if q <= 0.0 {
            broke = idx;
            break;
        }
        pqtable.insert(v.to_bits(), q as f32);
        pre_q = q;
        k += l;
    }
    for &v in &uniq[broke..] {
        pqtable.insert(v.to_bits(), 0.0);
    }
    pqtable
}

/// Find the summit (lower-median midpoint of the max-treatment segments) and,
/// if the peak is long enough and clears the cutoff, push it.
fn close_peak(
    tid: i32,
    content: &[(i32, i32, f32, f32, f32)], // (start, end, treat, ctrl, qscore)
    min_len: i32,
    cutoff: f32,
    cache: &mut HashMap<(i32, u32), f32>,
    pqtable: &HashMap<u32, f32>,
    peaks: &mut Vec<Peak>,
) {
    if content.is_empty() || content.last().unwrap().1 - content[0].0 < min_len {
        return;
    }
    let mut summit_value = 0f32;
    let mut tsummit: Vec<usize> = Vec::new();
    for (i, &(_, _, tp, _, _)) in content.iter().enumerate() {
        if summit_value == 0.0 || summit_value < tp {
            tsummit = vec![i];
            summit_value = tp;
        } else if (tp - summit_value).abs() < f32::EPSILON {
            tsummit.push(i);
        }
    }
    let idx = tsummit[tsummit.len().div_ceil(2) - 1];
    let (s, e, st, sc, sq) = content[idx];
    if cutoff > sq {
        return;
    }
    let summit = (s + e) / 2;
    let ps = pscore(st, sc, cache);
    let q = *pqtable.get(&ps.to_bits()).unwrap_or(&0.0);
    peaks.push(Peak {
        tid,
        start: content[0].0,
        end: content.last().unwrap().1,
        summit,
        pileup: st,
        pscore: ps,
        qscore: q,
        // MACS computes fold change in f64 (pseudocount is a C double) then
        // narrows to f32; doing it in f32 differs by 1 ulp at %.6g boundaries.
        fc: ((f64::from(st) + 1.0) / (f64::from(sc) + 1.0)) as f32,
    });
}

/// Per-chromosome peak calling over the combined `(end, treat, ctrl)` segments.
#[allow(clippy::too_many_arguments)]
fn call_chrom(
    tid: i32,
    segs: &[Seg],
    max_gap: i32,
    min_len: i32,
    cutoff: f32,
    cache: &mut HashMap<(i32, u32), f32>,
    pqtable: &HashMap<u32, f32>,
    peaks: &mut Vec<Peak>,
) {
    let mut content: Vec<(i32, i32, f32, f32, f32)> = Vec::new();
    let mut pre = 0i32;
    let mut lastp = 0i32;
    for &(end, tv, cv) in segs {
        let ps = pscore(tv, cv, cache);
        let q = *pqtable.get(&ps.to_bits()).unwrap_or(&0.0);
        let start = pre;
        pre = end;
        if q > cutoff {
            if !content.is_empty() && start - lastp > max_gap {
                close_peak(tid, &content, min_len, cutoff, cache, pqtable, peaks);
                content.clear();
            }
            content.push((start, end, tv, cv, q));
            lastp = end;
        }
    }
    close_peak(tid, &content, min_len, cutoff, cache, pqtable, peaks);
}

/// Run the no-control callpeak final stage: build the q-value table, call peaks
/// (qscore > -log10(0.05), gaps <= tsize, length >= d), keep fc >= 1.0.
#[must_use]
pub fn call_peaks(tags: &Tags, d: i32, treat: &[RawTrack], control: &[RawTrack]) -> Vec<Peak> {
    let max_gap = tags.tsize as i32;
    let min_len = d;
    let cutoff = (-(0.05f64).log10()) as f32;
    let ctrl_by_tid: HashMap<i32, (&[i32], &[f32])> = control
        .iter()
        .map(|(t, p, v)| (*t, (p.as_slice(), v.as_slice())))
        .collect();

    let combined: Vec<(i32, Vec<Seg>)> = treat
        .iter()
        .filter_map(|(tid, tp, tv)| {
            ctrl_by_tid
                .get(tid)
                .map(|(cp, cv)| (*tid, combine(tp, tv, cp, cv)))
        })
        .collect();

    let mut cache: HashMap<(i32, u32), f32> = HashMap::new();
    let all: Vec<&[Seg]> = combined.iter().map(|(_, c)| c.as_slice()).collect();
    let pqtable = build_pqtable(&all, &mut cache);

    let mut peaks = Vec::new();
    for (tid, segs) in &combined {
        call_chrom(
            *tid, segs, max_gap, min_len, cutoff, &mut cache, &pqtable, &mut peaks,
        );
    }
    peaks.retain(|p| p.fc >= 1.0);
    peaks
}

/// Format like Python `%.6g` (6 significant figures) for narrowPeak values.
fn py_g(v: f32) -> String {
    let v = f64::from(v);
    if v == 0.0 {
        return "0".to_string();
    }
    // Python decides fixed vs scientific from the exponent *after* rounding to 6
    // figures, so branch on the `{:.5e}` exponent rather than `floor(log10)` of
    // the raw value (they differ at the 1e-4 / 1e6 boundaries).
    let sci = format!("{v:.5e}");
    let (mantissa, e) = sci.split_once('e').unwrap();
    let exp: i32 = e.parse().unwrap();
    if (-4..6).contains(&exp) {
        let decimals = (5 - exp).max(0) as usize;
        let s = format!("{v:.decimals$}");
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s
        }
    } else {
        let mantissa = if mantissa.contains('.') {
            mantissa.trim_end_matches('0').trim_end_matches('.')
        } else {
            mantissa
        };
        format!(
            "{mantissa}e{}{:02}",
            if exp < 0 { '-' } else { '+' },
            exp.abs()
        )
    }
}

/// Write `<name>_peaks.narrowPeak` and `<name>_summits.bed` into `outdir`,
/// sorted by chromosome name then start and globally numbered from 1.
pub fn write_outputs(
    peaks: &[Peak],
    names: &HashMap<i32, String>,
    name: &str,
    outdir: &Path,
) -> Result<()> {
    let mut order: Vec<&Peak> = peaks.iter().collect();
    order.sort_by(|a, b| {
        let na = names.get(&a.tid).map_or("", String::as_str);
        let nb = names.get(&b.tid).map_or("", String::as_str);
        na.cmp(nb).then(a.start.cmp(&b.start))
    });
    let mut np = BufWriter::new(
        File::create(outdir.join(format!("{name}_peaks.narrowPeak"))).map_err(RsomicsError::Io)?,
    );
    let mut sm = BufWriter::new(
        File::create(outdir.join(format!("{name}_summits.bed"))).map_err(RsomicsError::Io)?,
    );
    for (i, p) in order.iter().enumerate() {
        let n = i + 1;
        let chrom = names.get(&p.tid).map_or("?", String::as_str);
        let score = (10.0 * p.qscore) as i32;
        writeln!(
            np,
            "{chrom}\t{}\t{}\t{name}_peak_{n}\t{score}\t.\t{}\t{}\t{}\t{}",
            p.start,
            p.end,
            py_g(p.fc),
            py_g(p.pscore),
            py_g(p.qscore),
            p.summit - p.start
        )?;
        writeln!(
            sm,
            "{chrom}\t{}\t{}\t{name}_peak_{n}\t{}",
            p.summit,
            p.summit + 1,
            py_g(p.qscore)
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::py_g;

    #[test]
    fn py_g_matches_python_6g() {
        assert_eq!(py_g(0.0), "0");
        assert_eq!(py_g(18.5), "18.5");
        assert_eq!(py_g(100.0), "100");
        assert_eq!(py_g(123.456), "123.456");
        assert_eq!(py_g(2_000_000.0), "2e+06");
        assert_eq!(py_g(2f32.powi(-20)), "9.53674e-07");
        // boundaries decided after rounding, not on floor(log10) of the raw value
        assert_eq!(py_g(999_999.5), "1e+06");
        assert_eq!(py_g(1e-4), "0.0001");
    }
}
