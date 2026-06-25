use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, Tool, ToolMeta};

use rsomics_macs::{CallPeakOpts, run_callpeak};

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

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    fn run_inner(self) -> Result<()> {
        run_callpeak(&CallPeakOpts {
            treatment: self.treatment,
            control: self.control,
            keep_dup: self.keep_dup,
            name: self.name,
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
