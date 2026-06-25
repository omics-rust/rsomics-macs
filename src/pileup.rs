//! Fragment-extension pileup and the no-control dynamic local lambda.
//!
//! Value-exact port of MACS3 `se_all_in_one_pileup` + `__call_peaks_wo_control`
//! per `.autopilot/state/macs-pileup-spec.md`.

use crate::tags::Tags;

/// One bedGraph segment `[start, end)` carrying a value.
pub type BedGraph = Vec<(i32, i32, f32)>;

/// `se_all_in_one_pileup`: + tags cover `[p-five, p+three)`, - tags cover
/// `[m-three, m+five)`; sweep emits the pileup *before* each event, floored at
/// `baseline` and scaled by `scale`.
fn se_pileup(
    plus: &[i32],
    minus: &[i32],
    five: i32,
    three: i32,
    scale: f32,
    baseline: f32,
    rlength: i32,
) -> (Vec<i32>, Vec<f32>) {
    let mut starts: Vec<i32> = Vec::with_capacity(plus.len() + minus.len());
    let mut ends: Vec<i32> = Vec::with_capacity(plus.len() + minus.len());
    for &p in plus {
        starts.push(p - five);
        ends.push(p + three);
    }
    for &m in minus {
        starts.push(m - three);
        ends.push(m + five);
    }
    starts.sort_unstable();
    ends.sort_unstable();
    clip(&mut starts, rlength);
    clip(&mut ends, rlength);

    let lx = starts.len();
    let mut rp = Vec::new();
    let mut rv = Vec::new();
    if lx == 0 {
        return (rp, rv);
    }
    let mut pileup = 0i32;
    let (mut is, mut ie) = (0usize, 0usize);
    let mut pre_p = starts[0].min(ends[0]);
    if pre_p != 0 {
        rp.push(pre_p);
        rv.push(0.0f32.max(baseline));
    }
    while is < lx && ie < lx {
        let s = starts[is];
        let e = ends[ie];
        if s < e {
            if s != pre_p {
                rp.push(s);
                rv.push((pileup as f32 * scale).max(baseline));
                pre_p = s;
            }
            pileup += 1;
            is += 1;
        } else if s > e {
            if e != pre_p {
                rp.push(e);
                rv.push((pileup as f32 * scale).max(baseline));
                pre_p = e;
            }
            pileup -= 1;
            ie += 1;
        } else {
            is += 1;
            ie += 1;
        }
    }
    while ie < lx {
        let e = ends[ie];
        if e != pre_p {
            rp.push(e);
            rv.push((pileup as f32 * scale).max(baseline));
            pre_p = e;
        }
        pileup -= 1;
        ie += 1;
    }
    (rp, rv)
}

/// MACS `fix_coordinates`: clamp leading negatives to 0 and trailing values
/// past the reference end to `rlength` (the array is already sorted).
fn clip(arr: &mut [i32], rlength: i32) {
    for v in arr.iter_mut() {
        if *v < 0 {
            *v = 0;
        } else {
            break;
        }
    }
    for v in arr.iter_mut().rev() {
        if *v > rlength {
            *v = rlength;
        } else {
            break;
        }
    }
}

/// Collapse the `(end_pos, value)` stream into bedGraph segments, merging runs
/// whose value differs by `<= 1e-5` (the MACS3 bedGraph writer's tolerance).
fn to_bedgraph(p: &[i32], v: &[f32]) -> BedGraph {
    let mut out = BedGraph::new();
    if p.is_empty() {
        return out;
    }
    let mut run_start = 0i32;
    let mut run_v = v[0];
    let mut prev_end = 0i32;
    for i in 0..p.len() {
        if (v[i] - run_v).abs() > 1e-5 {
            out.push((run_start, prev_end, run_v));
            run_start = prev_end;
            run_v = v[i];
        }
        prev_end = p[i];
    }
    out.push((run_start, prev_end, run_v));
    out
}

/// Truncate the treat/control bedGraph pair to their common right extent. MACS3
/// writes both tracks from one combined `over_two_pv_array` walk, which stops
/// when either array is exhausted — so the longer track's trailing segments are
/// dropped from the written bedGraph.
pub fn truncate_to_common(treat: &mut BedGraph, control: &mut BedGraph) {
    let m = treat
        .last()
        .map_or(0, |s| s.1)
        .min(control.last().map_or(0, |s| s.1));
    for bg in [treat, control] {
        bg.retain(|s| s.0 < m);
        if let Some(last) = bg.last_mut()
            && last.1 > m
        {
            last.1 = m;
        }
    }
}

fn sorted_tids(tags: &Tags) -> Vec<i32> {
    let mut tids: Vec<i32> = tags.plus.keys().chain(tags.minus.keys()).copied().collect();
    tids.sort_unstable();
    tids.dedup();
    tids
}

/// Treatment pileup bedGraph per reference: each tag extended to fragment
/// length `d` (+ tag `[p, p+d)`, - tag `[m-d, m)`).
#[must_use]
pub fn treat_pileup(tags: &Tags, d: i32) -> Vec<(i32, BedGraph)> {
    let empty = Vec::new();
    sorted_tids(tags)
        .iter()
        .map(|&tid| {
            let pt = tags.plus.get(&tid).unwrap_or(&empty);
            let mt = tags.minus.get(&tid).unwrap_or(&empty);
            let rlen = *tags.lengths.get(&tid).unwrap_or(&i32::MAX);
            let (p, v) = se_pileup(pt, mt, 0, d, 1.0, 0.0, rlen);
            (tid, to_bedgraph(&p, &v))
        })
        .collect()
}

/// No-control dynamic lambda bedGraph per reference: the treatment tags piled
/// with a symmetric `llocal` window scaled by `d/llocal`, floored at
/// `lambda_bg = d * total / gsize`.
#[must_use]
pub fn control_lambda(tags: &Tags, d: i32, gsize: f64, llocal: i32) -> Vec<(i32, BedGraph)> {
    let lambda_bg = (f64::from(d) * tags.total as f64 / gsize) as f32;
    let scale = d as f32 / llocal as f32;
    let half = llocal / 2;
    let (five, three) = (half, llocal - half);
    let empty = Vec::new();
    sorted_tids(tags)
        .iter()
        .map(|&tid| {
            let pt = tags.plus.get(&tid).unwrap_or(&empty);
            let mt = tags.minus.get(&tid).unwrap_or(&empty);
            let rlen = *tags.lengths.get(&tid).unwrap_or(&i32::MAX);
            let (p, v) = se_pileup(pt, mt, five, three, scale, lambda_bg, rlen);
            (tid, to_bedgraph(&p, &v))
        })
        .collect()
}
