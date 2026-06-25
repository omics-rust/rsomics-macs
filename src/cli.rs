use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};

use rsomics_macs::{CallPeakOpts, run_callpeak};

/// Effective genome size from a number or an `hs`/`mm`/`ce`/`dm` shortcut.
fn parse_gsize(s: &str) -> Result<f64> {
    Ok(match s {
        "hs" => 2_913_022_398.0,
        "mm" => 2_652_783_500.0,
        "ce" => 100_286_401.0,
        "dm" => 142_573_017.0,
        _ => s
            .parse::<f64>()
            .map_err(|_| RsomicsError::InvalidInput(format!("invalid --gsize: {s}")))?,
    })
}

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(
    name = "rsomics-macs",
    version,
    about = "Model-based ChIP-seq peak caller — MACS3 callpeak port (work in progress)",
    long_about = None
)]
pub struct Cli {
    /// ChIP-seq treatment BAM file(s). REQUIRED. (`-t` is reserved by the
    /// shared `--threads` flag; treatment is long-only for now.)
    #[arg(long = "treatment", required = true, num_args = 1..)]
    pub treatment: Vec<PathBuf>,

    /// Control (input) BAM file(s).
    #[arg(short = 'c', long = "control", num_args = 0..)]
    pub control: Vec<PathBuf>,

    /// Duplicate handling: an integer cap per position, "all", or "auto".
    #[arg(long = "keep-dup", default_value = "1")]
    pub keep_dup: String,

    /// Experiment name, used as the output file prefix.
    #[arg(short = 'n', long = "name", default_value = "NA")]
    pub name: String,

    /// Effective genome size: a number (e.g. 1.0e9) or shortcut hs/mm/ce/dm.
    #[arg(short = 'g', long = "gsize", default_value = "hs")]
    pub gsize: String,

    /// Skip the shifting-model step and extend reads to `--extsize`.
    #[arg(long = "nomodel")]
    pub nomodel: bool,

    /// Fragment size used with `--nomodel`.
    #[arg(long = "extsize", default_value_t = 200)]
    pub extsize: i32,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    fn run_inner(self) -> Result<()> {
        let gsize = parse_gsize(&self.gsize)?;
        run_callpeak(&CallPeakOpts {
            treatment: self.treatment,
            control: self.control,
            keep_dup: self.keep_dup,
            name: self.name,
            gsize,
            nomodel: self.nomodel,
            extsize: self.extsize,
        })
    }
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        self.run_inner()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
