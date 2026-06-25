//! MACS3 `PeakModel.build()` — fragment-size `d` estimation.
//!
//! Naive pileup-based peak finding on each strand, +/- peak pairing, then a
//! cross-correlation of the paired-peak tag profiles whose dominant lag is the
//! estimated fragment length `d`.

use rsomics_common::{Result, RsomicsError};

use crate::tags::Tags;

/// Result of model building.
pub struct PeakModel {
    pub d: i32,
    pub alternative_d: Vec<i32>,
    pub scan_window: i32,
}

/// Python `round()` — round half to even — then truncate to integer.
fn banker_round(x: f64) -> i64 {
    let f = x.floor();
    let diff = x - f;
    if diff > 0.5 {
        f as i64 + 1
    } else if diff < 0.5 {
        f as i64
    } else {
        let fi = f as i64;
        if fi % 2 == 0 { fi } else { fi + 1 }
    }
}

/// `naive_quick_pileup`: extend each sorted 5' position by `±ext` and sweep into
/// a bedGraph `(end_pos, value)` stream (value = integer pileup as f32).
fn naive_quick_pileup(poss: &[i32], ext: i32) -> (Vec<i32>, Vec<f32>) {
    let n = poss.len();
    let mut rp = Vec::new();
    let mut rv = Vec::new();
    if n == 0 {
        return (rp, rv);
    }
    let starts: Vec<i32> = poss.iter().map(|&p| (p - ext).max(0)).collect();
    let ends: Vec<i32> = poss.iter().map(|&p| p + ext).collect();
    let mut is = 0usize;
    let mut ie = 0usize;
    let mut pileup = 0i32;
    let mut pre_p = starts[0].min(ends[0]);
    if pre_p != 0 {
        rp.push(pre_p);
        rv.push(0.0);
    }
    while is < n && ie < n {
        let s = starts[is];
        let e = ends[ie];
        if s < e {
            if s != pre_p {
                rp.push(s);
                rv.push(pileup as f32);
            }
            pileup += 1;
            pre_p = s;
            is += 1;
        } else if e < s {
            if e != pre_p {
                rp.push(e);
                rv.push(pileup as f32);
            }
            pileup -= 1;
            pre_p = e;
            ie += 1;
        } else {
            is += 1;
            ie += 1;
        }
    }
    while ie < n {
        let e = ends[ie];
        if e != pre_p {
            rp.push(e);
            rv.push(pileup as f32);
        }
        pileup -= 1;
        pre_p = e;
        ie += 1;
    }
    (rp, rv)
}

/// `__close_peak`: pick the lower-median midpoint of the maximum-height segments
/// and keep the peak iff its summit height is `< max_v`.
fn close_peak(content: &[(i32, i32, f32)], peaks: &mut Vec<(i32, f32)>, max_v: f64, min_len: i32) {
    if content.is_empty() {
        return;
    }
    let peak_len = content.last().unwrap().1 - content[0].0;
    if peak_len < min_len {
        return;
    }
    let mut summit_value = 0f64;
    let mut tsummit: Vec<i32> = Vec::new();
    for &(ts, te, tv) in content {
        let tvf = f64::from(tv);
        if summit_value == 0.0 || summit_value < tvf {
            tsummit = vec![(ts + te) / 2];
            summit_value = tvf;
        } else if (tvf - summit_value).abs() < f64::EPSILON {
            tsummit.push((ts + te) / 2);
        }
    }
    if summit_value < max_v {
        let n = tsummit.len();
        let idx = ((n + 1) as f64 / 2.0) as usize - 1;
        peaks.push((tsummit[idx], summit_value as f32));
    }
}

/// `naive_call_peaks`: collect above-`min_v` segments into peaks, merging within
/// `max_gap`, keeping those at least `min_len` long.
fn naive_call_peaks(
    p: &[i32],
    v: &[f32],
    min_v: f64,
    max_v: f64,
    max_gap: i32,
    min_len: i32,
) -> Vec<(i32, f32)> {
    let mut peaks = Vec::new();
    let mut content: Vec<(i32, i32, f32)> = Vec::new();
    let mut pre_p = 0i32;
    for i in 0..p.len() {
        let cur_p = p[i];
        let cur_v = f64::from(v[i]);
        if cur_v > min_v {
            if content.is_empty() {
                content.push((pre_p, cur_p, v[i]));
            } else {
                let gap = pre_p - content.last().unwrap().1;
                if gap <= max_gap {
                    content.push((pre_p, cur_p, v[i]));
                } else {
                    close_peak(&content, &mut peaks, max_v, min_len);
                    content.clear();
                    content.push((pre_p, cur_p, v[i]));
                }
            }
        }
        pre_p = cur_p;
    }
    close_peak(&content, &mut peaks, max_v, min_len);
    peaks
}

fn naive_find_peaks(tags: &[i32], peaksize: i32, min_v: f64, max_v: f64) -> Vec<(i32, f32)> {
    if tags.is_empty() {
        return Vec::new();
    }
    let (p, v) = naive_quick_pileup(tags, peaksize / 2);
    naive_call_peaks(&p, &v, min_v, max_v, 50, 200)
}

/// `__find_pair_center`: pair a + peak with a downstream - peak of comparable
/// height; the paired centre is the midpoint.
fn find_pair_center(plus: &[(i32, f32)], minus: &[(i32, f32)], peaksize: i32) -> Vec<i32> {
    let mut centers = Vec::new();
    let (mut ip, mut im) = (0usize, 0usize);
    let mut im_prev = 0usize;
    let mut flag_overlap = false;
    while ip < plus.len() && im < minus.len() {
        let (pp, pn) = plus[ip];
        let (mp, mn) = minus[im];
        if pp - peaksize > mp {
            im += 1;
        } else if pp + peaksize < mp {
            ip += 1;
            im = im_prev;
            flag_overlap = false;
        } else {
            if !flag_overlap {
                flag_overlap = true;
                im_prev = im;
            }
            let r = f64::from(pn) / f64::from(mn);
            if r < 2.0 && r > 0.5 && pp < mp {
                centers.push((pp + mp) / 2);
            }
            im += 1;
        }
    }
    centers
}

/// `__model_add_line`: project tags within `±605` of each paired centre into the
/// `start`/`end` difference arrays.
fn add_line(centers: &[i32], tags: &[i32], start: &mut [i32], end: &mut [i32]) {
    const PSIZE_ADJ: i32 = 605;
    let mut im = 0usize;
    for &p1 in centers {
        while im < tags.len() && tags[im] < p1 - PSIZE_ADJ {
            im += 1;
        }
        let mut j = im;
        while j < tags.len() && tags[j] <= p1 + PSIZE_ADJ {
            let p2 = tags[j];
            let s = (p2 - p1 + 600).max(0) as usize;
            let e = (p2 - p1 + 610).min(1210) as usize;
            start[s] += 1;
            end[e] -= 1;
            j += 1;
        }
    }
}

fn count_line(start: &[i32], end: &[i32]) -> Vec<i32> {
    let mut line = vec![0i32; start.len()];
    let mut pileup = 0i64;
    for i in 0..start.len() {
        pileup += i64::from(start[i]) + i64::from(end[i]);
        line[i] = pileup as i32;
    }
    line
}

/// `(line - mean) / (std * len)` with population std (ddof=0).
fn normalize(line: &[i32]) -> Vec<f64> {
    let n = line.len() as f64;
    let mean = line.iter().map(|&x| f64::from(x)).sum::<f64>() / n;
    let var = line
        .iter()
        .map(|&x| {
            let d = f64::from(x) - mean;
            d * d
        })
        .sum::<f64>()
        / n;
    let denom = var.sqrt() * n;
    line.iter()
        .map(|&x| (f64::from(x) - mean) / denom)
        .collect()
}

/// `np.correlate(a, b, "full")[ws-peaksize : ws+peaksize]`.
fn correlate_slice(a: &[f64], b: &[f64], peaksize: usize) -> Vec<f64> {
    let n = a.len() as i64;
    let center = n - 1;
    let ws = a.len();
    let mut out = Vec::with_capacity(2 * peaksize);
    for j in (ws - peaksize)..(ws + peaksize) {
        let lag = j as i64 - center;
        let i_lo = lag.max(0);
        let i_hi = (n + lag).min(n);
        let mut s = 0f64;
        let mut i = i_lo;
        while i < i_hi {
            s += a[i as usize] * b[(i - lag) as usize];
            i += 1;
        }
        out.push(s);
    }
    out
}

/// `smooth(x, window="flat", window_len=11)`: reflected pad, 11-wide moving
/// average, trimmed back to the input length.
fn smooth(x: &[f64]) -> Vec<f64> {
    const WL: usize = 11;
    let n = x.len();
    let mut s = Vec::with_capacity(n + 2 * (WL - 1));
    for i in (1..WL).rev() {
        s.push(x[i]);
    }
    s.extend_from_slice(x);
    for i in 0..(WL - 1) {
        s.push(x[n - 1 - i]);
    }
    let mut y = Vec::with_capacity(s.len() - WL + 1);
    for i in 0..(s.len() - WL + 1) {
        let sum: f64 = s[i..i + WL].iter().sum();
        y.push(sum / WL as f64);
    }
    y[(WL / 2)..(y.len() - WL / 2)].to_vec()
}

pub fn build(
    tags: &Tags,
    gsize: f64,
    bw: i32,
    lmfold: f64,
    umfold: f64,
    d_min: i32,
) -> Result<PeakModel> {
    let peaksize = 2 * bw;
    let total = tags.total as f64;
    let min_tags = banker_round(total * lmfold * f64::from(peaksize) / gsize / 2.0) as f64;
    let max_tags = banker_round(total * umfold * f64::from(peaksize) / gsize / 2.0) as f64;

    let mut tids: Vec<i32> = tags.plus.keys().chain(tags.minus.keys()).copied().collect();
    tids.sort_unstable();
    tids.dedup();

    let empty: Vec<i32> = Vec::new();
    let mut paired: Vec<(i32, Vec<i32>)> = Vec::new();
    let mut total_pairs = 0usize;
    for &tid in &tids {
        let pt = tags.plus.get(&tid).unwrap_or(&empty);
        let mt = tags.minus.get(&tid).unwrap_or(&empty);
        let pp = naive_find_peaks(pt, peaksize, min_tags, max_tags);
        let mp = naive_find_peaks(mt, peaksize, min_tags, max_tags);
        let centers = find_pair_center(&pp, &mp, peaksize);
        total_pairs += centers.len();
        paired.push((tid, centers));
    }
    if total_pairs < 100 {
        return Err(RsomicsError::InvalidInput(format!(
            "only {total_pairs} paired peaks (need >=100); use --nomodel --extsize"
        )));
    }

    let ws = 1 + 2 * peaksize as usize + 10;
    let mut ps = vec![0i32; ws];
    let mut pe = vec![0i32; ws];
    let mut ms = vec![0i32; ws];
    let mut me = vec![0i32; ws];
    for (tid, centers) in &paired {
        let pt = tags.plus.get(tid).unwrap_or(&empty);
        let mt = tags.minus.get(tid).unwrap_or(&empty);
        add_line(centers, pt, &mut ps, &mut pe);
        add_line(centers, mt, &mut ms, &mut me);
    }
    let plus_data = normalize(&count_line(&ps, &pe));
    let minus_data = normalize(&count_line(&ms, &me));

    let ycorr = smooth(&correlate_slice(&minus_data, &plus_data, peaksize as usize));
    let nlen = ycorr.len();
    let half = (nlen / 2) as f64;
    let xcorr: Vec<f64> = (0..nlen)
        .map(|i| -half + i as f64 * (2.0 * half) / ((nlen - 1) as f64))
        .collect();

    let mut maxima: Vec<usize> = (1..nlen - 1)
        .filter(|&i| {
            ycorr[i - 1] < ycorr[i] && ycorr[i] > ycorr[i + 1] && xcorr[i] > f64::from(d_min)
        })
        .collect();
    if maxima.is_empty() {
        return Err(RsomicsError::InvalidInput(
            "model cross-correlation has no local maximum past d_min".into(),
        ));
    }
    maxima.sort_by(|&a, &b| ycorr[b].partial_cmp(&ycorr[a]).unwrap());
    let d = xcorr[maxima[0]] as i32;
    let mut alternative_d: Vec<i32> = maxima.iter().map(|&i| xcorr[i] as i32).collect();
    alternative_d.sort_unstable();
    let scan_window = d.max(10) * 2;
    Ok(PeakModel {
        d,
        alternative_d,
        scan_window,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn banker_round_half_to_even() {
        assert_eq!(banker_round(24.0), 24);
        assert_eq!(banker_round(2.5), 2);
        assert_eq!(banker_round(3.5), 4);
        assert_eq!(banker_round(24.4), 24);
        assert_eq!(banker_round(24.6), 25);
    }
}
