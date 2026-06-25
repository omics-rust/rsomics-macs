//! Phase C differential: treatment pileup + no-control dynamic lambda must match
//! MACS3 3.0.4 `callpeak --bdg --nomodel --extsize 144 -g 1000000` on the
//! committed ChIP fixture. Always-run (no macs3 needed); the goldens are macs3
//! 3.0.4's `*_treat_pileup.bdg` / `*_control_lambda.bdg`.

use std::path::Path;

use rsomics_macs::{pileup, tags};

const GOLDEN: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/golden");

fn parse_bdg(text: &str) -> Vec<(i32, i32, f32)> {
    text.lines()
        .filter_map(|l| {
            let mut it = l.split('\t');
            let _chrom = it.next()?;
            let s: i32 = it.next()?.parse().ok()?;
            let e: i32 = it.next()?.parse().ok()?;
            let v: f32 = it.next()?.parse().ok()?;
            Some((s, e, v))
        })
        .collect()
}

fn assert_bdg_eq(ours: &[(i32, i32, f32)], theirs: &[(i32, i32, f32)], what: &str) {
    assert_eq!(
        ours.len(),
        theirs.len(),
        "{what}: segment count ours={} macs3={}",
        ours.len(),
        theirs.len()
    );
    for (i, (o, t)) in ours.iter().zip(theirs).enumerate() {
        assert_eq!(o.0, t.0, "{what} seg {i}: start ours={} macs3={}", o.0, t.0);
        assert_eq!(o.1, t.1, "{what} seg {i}: end ours={} macs3={}", o.1, t.1);
        assert!(
            (o.2 - t.2).abs() <= 1e-4,
            "{what} seg {i}: value ours={} macs3={}",
            o.2,
            t.2
        );
    }
}

fn loaded() -> tags::Tags {
    let bam = Path::new(GOLDEN).join("chip_sim.bam");
    let mut t = tags::load_bam(&bam).expect("load chip_sim.bam");
    t.filter_dup(1);
    t
}

fn golden(name: &str) -> Vec<(i32, i32, f32)> {
    parse_bdg(&std::fs::read_to_string(Path::new(GOLDEN).join(name)).expect("read golden"))
}

#[test]
fn bedgraphs_match_macs3() {
    let t = loaded();
    let mut treat = pileup::treat_pileup(&t, 144).remove(0).1;
    let mut control = pileup::control_lambda(&t, 144, 1_000_000.0, 10000)
        .remove(0)
        .1;
    // MACS3 writes both tracks from one combined walk that drops the longer
    // track's trailing segments.
    pileup::truncate_to_common(&mut treat, &mut control);
    assert_bdg_eq(&treat, &golden("chip_sim.treat_pileup.bdg"), "treat_pileup");
    assert_bdg_eq(
        &control,
        &golden("chip_sim.control_lambda.bdg"),
        "control_lambda",
    );
}
