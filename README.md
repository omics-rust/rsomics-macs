# rsomics-macs

Model-based ChIP-seq peak caller — a Rust port of MACS3 `callpeak` for the
single-end, no-control case. Output (`narrowPeak` + `summits.bed`) is
**byte-identical** to MACS3 3.0.4, at roughly **10× the throughput** and **a
third of the memory**.

```sh
cargo install rsomics-macs
```

## Usage

```sh
# predict the fragment length, then call peaks
rsomics-macs --treatment chip.bam -g hs -n mysample --outdir peaks/

# skip modelling and extend reads to a fixed fragment size
rsomics-macs --treatment chip.bam -g 1.0e9 --nomodel --extsize 200 -n mysample
```

Writes `<name>_peaks.narrowPeak` and `<name>_summits.bed` into `--outdir`.

`-g` takes an effective genome size as a number or one of the shortcuts `hs`,
`mm`, `ce`, `dm`. `--keep-dup` takes an integer cap per position (default `1`) or
`all`.

## Scope

This release does single-end, no-control peak calling — the dynamic local Poisson
background is built from the treatment's 10 kb window (MACS's no-control mode,
where the small local window is omitted because the ChIP signal would inflate it).

Passing `-c`/`--control` errors loudly rather than silently ignoring it; control-
based calling (the three-window `d`/slocal/llocal background) is a later increment.

## Performance

On a 400 000-read / 2000-peak fixture (mini_m2, Apple Silicon, single-threaded),
versus `macs3 3.0.4`:

| | ours | macs3 3.0.4 | ratio |
|---|---|---|---|
| CPU (`--nomodel`) | 146 ms | 1.44 s | 9.86× faster |
| CPU (model path) | ~170 ms | ~1.68 s | 9.61× faster |
| peak RSS | 43 MB | 124 MB | 2.87× smaller |

`narrowPeak` and `summits.bed` are byte-identical to macs3 on every fixture
tested, including the model path (predicted `d` matches exactly). The end-to-end
byte comparison runs in CI against a checked-in macs3 golden (`tests/compat.rs`).

## Origin

This crate is an independent Rust reimplementation of MACS3 `callpeak`, informed
by:

- The published method: Zhang et al., *Model-based Analysis of ChIP-Seq (MACS)*,
  Genome Biology 2008, [doi:10.1186/gb-2008-9-9-r137](https://doi.org/10.1186/gb-2008-9-9-r137).
- The MACS3 source (BSD-3-Clause, a permissive license that allows reading and
  citing): the SAM/BAM tag parser, `se_all_in_one_pileup`, the no-control dynamic
  lambda, `log10_poisson_cdf_Q`, the Benjamini–Hochberg q-table, and the
  peak-closing / summit-selection logic.
- Black-box behaviour testing against the `macs3` binary.

Test fixtures are independently generated (`tests/golden/simulate_chip.py`).

License: MIT OR Apache-2.0.
Upstream credit: [MACS3](https://github.com/macs3-project/MACS) (BSD-3-Clause).
