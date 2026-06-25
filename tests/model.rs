//! Phase B differential: the PeakModel fragment-length `d` must match MACS3
//! 3.0.4 `callpeak` on the committed synthetic ChIP fixture. Always-run (no
//! macs3 needed) — the golden values are macs3 3.0.4's reported numbers.

use std::path::Path;

use rsomics_macs::{model, tags};

const GOLDEN: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/golden");

#[test]
fn model_d_matches_macs3() {
    let bam = Path::new(GOLDEN).join("chip_sim.bam");
    let mut t = tags::load_bam(&bam).expect("load chip_sim.bam");
    assert_eq!(t.tsize, 50, "tag size");
    t.filter_dup(1);
    // macs3 callpeak --keep-dup 1: tags after filtering = 15970.
    assert_eq!(t.total, 15970, "tags after dup-filter");
    let m = model::build(&t, 1_000_000.0, 300, 5.0, 50.0, 20).expect("model builds");
    // macs3 callpeak -g 1000000: predicted fragment length 144, alternatives 144.
    assert_eq!(m.d, 144, "predicted fragment length d");
    assert_eq!(m.alternative_d, vec![144], "alternative d");
}
