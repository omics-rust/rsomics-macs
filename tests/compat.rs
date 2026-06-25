//! End-to-end compat: the binary's `narrowPeak` and `summits.bed` must be
//! byte-identical to MACS3 3.0.4 `callpeak --nomodel --extsize 144 -g 1000000`
//! on the committed ChIP fixture. Always-run — the goldens are checked-in macs3
//! output, so no macs3 install is needed.

use std::path::Path;
use std::process::Command;

const GOLDEN: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/golden");

#[test]
fn binary_output_byte_identical_to_macs3() {
    let bam = Path::new(GOLDEN).join("chip_sim.bam");
    let out = tempfile::tempdir().expect("tempdir");
    let status = Command::new(env!("CARGO_BIN_EXE_rsomics-macs"))
        .arg("--treatment")
        .arg(&bam)
        .args(["-g", "1000000", "--nomodel", "--extsize", "144", "-n", "p"])
        .arg("--outdir")
        .arg(out.path())
        .status()
        .expect("run rsomics-macs");
    assert!(status.success(), "binary exited non-zero");

    for (produced, golden) in [
        ("p_peaks.narrowPeak", "chip_sim.narrowPeak"),
        ("p_summits.bed", "chip_sim.summits.bed"),
    ] {
        let got = std::fs::read_to_string(out.path().join(produced)).expect(produced);
        let want = std::fs::read_to_string(Path::new(GOLDEN).join(golden)).expect(golden);
        assert_eq!(got, want, "{produced} diverges from macs3 golden {golden}");
    }
}
