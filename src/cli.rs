use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_sc_cell_cycle::{CellCycleParams, open_output, run};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-sc-cell-cycle", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// 10x MTX directory (matrix.mtx[.gz] + features.tsv[.gz] + barcodes.tsv[.gz]).
    pub input: PathBuf,

    #[arg(long = "s-genes")]
    s_genes: PathBuf,

    #[arg(long = "g2m-genes")]
    g2m_genes: PathBuf,

    #[arg(short = 'o', long, default_value = "-")]
    output: String,

    #[arg(long = "n-bins", default_value_t = 25)]
    n_bins: usize,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        self.common.install_rayon_pool()?;
        let params = CellCycleParams {
            n_bins: self.n_bins,
            seed: self.common.seed.unwrap_or(0) as u32,
        };
        let out = open_output(&self.output)?;
        let cells = run(&self.input, &self.s_genes, &self.g2m_genes, &params, out)?;
        if !self.common.quiet {
            eprintln!("scored {cells} cells");
        }
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Per-cell cell-cycle phase: S_score, G2M_score, and a phase call (S / G2M / G1).",
    origin: Some(Origin {
        upstream: "scanpy sc.tl.score_genes_cell_cycle",
        upstream_license: "BSD-3-Clause",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1038/nbt.3192"),
    }),
    usage_lines: &["<10x-mtx-dir> --s-genes <s.txt> --g2m-genes <g2m.txt> [-o out.tsv]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: None,
                long: "s-genes",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: true,
                default: None,
                description: "S-phase gene-list file, one gene id (features.tsv column 1) per line.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "g2m-genes",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: true,
                default: None,
                description: "G2M-phase gene-list file, one gene id per line.",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("String"),
                required: false,
                default: Some("-"),
                description: "Output TSV (cell_id<TAB>S_score<TAB>G2M_score<TAB>phase); '-' for stdout.",
                why_default: Some("Streams to stdout for pipeline composition."),
            },
            FlagSpec {
                short: None,
                long: "n-bins",
                aliases: &[],
                value: Some("<n>"),
                type_hint: Some("usize"),
                required: false,
                default: Some("25"),
                description: "Number of expression-level bins for control sampling.",
                why_default: Some("scanpy default n_bins=25."),
            },
            FlagSpec {
                short: None,
                long: "seed",
                aliases: &[],
                value: Some("<n>"),
                type_hint: Some("u64"),
                required: false,
                default: Some("0"),
                description: "Seed for the control-gene draw (numpy RandomState).",
                why_default: Some("scanpy default random_state=0; deterministic."),
            },
        ],
    }],
    examples: &[Example {
        description: "assign cell-cycle phase from Tirosh S/G2M gene lists",
        command: "rsomics-sc-cell-cycle filtered_feature_bc_matrix/ --s-genes s.txt --g2m-genes g2m.txt -o phase.tsv",
    }],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
