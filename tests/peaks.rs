//! Phase D/E differential: the called peaks (narrowPeak fields) must match MACS3
//! 3.0.4 `callpeak --nomodel --extsize 144 -g 1000000 -q 0.05` on the committed
//! ChIP fixture. Always-run; golden = macs3 3.0.4's narrowPeak (200 peaks).

use std::path::Path;

use rsomics_macs::{peaks, pileup, tags};

const GOLDEN: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/golden");

/// narrowPeak row: (start, end, score, fc, pscore, qscore, summit_offset).
struct Np {
    start: i32,
    end: i32,
    score: i64,
    fc: f64,
    pscore: f64,
    qscore: f64,
    offset: i32,
}

fn parse_narrowpeak(text: &str) -> Vec<Np> {
    text.lines()
        .filter_map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            if c.len() < 10 {
                return None;
            }
            Some(Np {
                start: c[1].parse().ok()?,
                end: c[2].parse().ok()?,
                score: c[4].parse().ok()?,
                fc: c[6].parse().ok()?,
                pscore: c[7].parse().ok()?,
                qscore: c[8].parse().ok()?,
                offset: c[9].parse().ok()?,
            })
        })
        .collect()
}

#[test]
fn narrowpeak_matches_macs3() {
    let bam = Path::new(GOLDEN).join("chip_sim.bam");
    let mut t = tags::load_bam(&bam).expect("load");
    t.filter_dup(1);
    let treat = pileup::treat_pileup_raw(&t, 144);
    let control = pileup::control_lambda_raw(&t, 144, 1_000_000.0, 10000);
    let ours = peaks::call_peaks(&t, 144, &treat, &control);

    let golden = parse_narrowpeak(
        &std::fs::read_to_string(Path::new(GOLDEN).join("chip_sim.narrowPeak")).expect("golden"),
    );
    assert_eq!(ours.len(), golden.len(), "peak count");
    for (i, (o, g)) in ours.iter().zip(&golden).enumerate() {
        assert_eq!(o.start, g.start, "peak {i} start");
        assert_eq!(o.end, g.end, "peak {i} end");
        assert_eq!(o.summit - o.start, g.offset, "peak {i} summit offset");
        assert_eq!(
            i64::from((10.0 * o.qscore) as i32),
            g.score,
            "peak {i} score=int(10*q)"
        );
        assert!(
            (f64::from(o.fc) - g.fc).abs() <= 1e-3,
            "peak {i} fc ours={} macs3={}",
            o.fc,
            g.fc
        );
        assert!(
            (f64::from(o.pscore) - g.pscore).abs() <= 1e-3,
            "peak {i} pscore ours={} macs3={}",
            o.pscore,
            g.pscore
        );
        assert!(
            (f64::from(o.qscore) - g.qscore).abs() <= 1e-3,
            "peak {i} qscore ours={} macs3={}",
            o.qscore,
            g.qscore
        );
    }
}
