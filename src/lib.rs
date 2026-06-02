use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use rsomics_common::{Result, RsomicsError};
use rsomics_sc_score_genes::{CountMatrix, ScoreParams, read_10x, read_gene_list, score};

pub struct CellCycle {
    pub barcodes: Vec<String>,
    pub s_score: Vec<f64>,
    pub g2m_score: Vec<f64>,
    pub phase: Vec<&'static str>,
}

pub struct CellCycleParams {
    pub n_bins: usize,
    pub seed: u32,
}

/// scanpy precedence: default S, then G2M where it outscores S, then G1 where
/// both scores are negative (the G1 rule is applied last and wins).
fn assign_phase(s: f64, g2m: f64) -> &'static str {
    if s < 0.0 && g2m < 0.0 {
        "G1"
    } else if g2m > s {
        "G2M"
    } else {
        "S"
    }
}

/// score_genes twice (ctrl_size = min over the two sets, as scanpy does), then
/// the deterministic phase call.
pub fn score_cell_cycle(
    m: &CountMatrix,
    s_genes: &[String],
    g2m_genes: &[String],
    params: &CellCycleParams,
) -> Result<CellCycle> {
    let ctrl_size = s_genes.len().min(g2m_genes.len());
    let mk = |names: &[String]| -> Result<Vec<f64>> {
        score(
            m,
            names,
            &ScoreParams {
                ctrl_size,
                n_bins: params.n_bins,
                seed: params.seed,
            },
        )
    };
    let s_score = mk(s_genes)?;
    let g2m_score = mk(g2m_genes)?;

    let phase = s_score
        .iter()
        .zip(&g2m_score)
        .map(|(&s, &g)| assign_phase(s, g))
        .collect();

    Ok(CellCycle {
        barcodes: m.barcodes.clone(),
        s_score,
        g2m_score,
        phase,
    })
}

pub fn write_tsv(cc: &CellCycle, out: impl Write) -> Result<()> {
    let mut w = BufWriter::with_capacity(1 << 20, out);
    w.write_all(b"cell_id\tS_score\tG2M_score\tphase\n")
        .map_err(RsomicsError::Io)?;
    let mut fmt = ryu::Buffer::new();
    let mut line: Vec<u8> = Vec::with_capacity(96);
    for i in 0..cc.barcodes.len() {
        line.clear();
        line.extend_from_slice(cc.barcodes[i].as_bytes());
        line.push(b'\t');
        line.extend_from_slice(fmt.format(cc.s_score[i]).as_bytes());
        line.push(b'\t');
        line.extend_from_slice(fmt.format(cc.g2m_score[i]).as_bytes());
        line.push(b'\t');
        line.extend_from_slice(cc.phase[i].as_bytes());
        line.push(b'\n');
        w.write_all(&line).map_err(RsomicsError::Io)?;
    }
    w.flush().map_err(RsomicsError::Io)?;
    Ok(())
}

pub fn open_output(path: &str) -> Result<Box<dyn Write>> {
    if path == "-" {
        Ok(Box::new(std::io::stdout().lock()))
    } else {
        Ok(Box::new(
            File::create(PathBuf::from(path)).map_err(RsomicsError::Io)?,
        ))
    }
}

pub fn run(
    mtx_dir: &Path,
    s_genes: &Path,
    g2m_genes: &Path,
    params: &CellCycleParams,
    out: impl Write,
) -> Result<usize> {
    let m = read_10x(mtx_dir)?;
    let s = read_gene_list(s_genes)?;
    let g2m = read_gene_list(g2m_genes)?;
    let cc = score_cell_cycle(&m, &s, &g2m, params)?;
    write_tsv(&cc, out)?;
    Ok(m.n_cells)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rsomics_sc_score_genes::{CountMatrix, Entry};

    fn matrix(n_genes: usize, n_cells: usize, triples: &[(u32, u32, f64)]) -> CountMatrix {
        CountMatrix {
            n_genes,
            n_cells,
            gene_ids: (0..n_genes).map(|g| format!("G{g}")).collect(),
            barcodes: (0..n_cells).map(|c| format!("C{c}")).collect(),
            entries: triples
                .iter()
                .map(|&(gene, cell, value)| Entry { gene, cell, value })
                .collect(),
        }
    }

    #[test]
    fn phase_precedence() {
        assert_eq!(assign_phase(0.5, -0.2), "S"); // S positive, G2M loses
        assert_eq!(assign_phase(0.1, 0.4), "G2M"); // G2M outscores S, both >= 0
        assert_eq!(assign_phase(-0.3, -0.1), "G1"); // both negative
        assert_eq!(assign_phase(0.0, 0.0), "S"); // tie at zero stays S
        assert_eq!(assign_phase(-0.5, 0.2), "G2M"); // only S<0, so not G1; G2M>S -> G2M
    }

    #[test]
    fn g2m_strictly_greater_for_g2m_call() {
        // when G2M does not strictly exceed S, phase stays S.
        assert_eq!(assign_phase(0.3, 0.3), "S");
        assert_eq!(assign_phase(0.3, 0.3000001), "G2M");
    }

    #[test]
    fn end_to_end_runs_and_calls_phases() {
        // 60 genes / 8 cells so control bins are populated after excluding the
        // 3+3 list genes; just checks shape + a deterministic phase per cell.
        let mut triples = Vec::new();
        for g in 0..60u32 {
            for c in 0..8u32 {
                let v = ((g * 7 + c * 3) % 11) as f64 + 1.0;
                triples.push((g, c, v));
            }
        }
        let m = matrix(60, 8, &triples);
        let s = vec!["G0".to_string(), "G1".to_string(), "G2".to_string()];
        let g2m = vec!["G30".to_string(), "G31".to_string(), "G32".to_string()];
        let params = CellCycleParams {
            n_bins: 25,
            seed: 0,
        };
        let cc = score_cell_cycle(&m, &s, &g2m, &params).unwrap();
        assert_eq!(cc.phase.len(), 8);
        for p in &cc.phase {
            assert!(matches!(*p, "S" | "G2M" | "G1"));
        }
    }
}
