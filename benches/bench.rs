use std::hint::black_box;
use std::path::PathBuf;

use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_sc_cell_cycle::{CellCycleParams, score_cell_cycle};
use rsomics_sc_score_genes::{read_10x, read_gene_list};

fn bench_cell_cycle(c: &mut Criterion) {
    let dir = match std::env::var("RSOMICS_BENCH_MTX") {
        Ok(d) => PathBuf::from(d),
        Err(_) => return,
    };
    let s = PathBuf::from(std::env::var("RSOMICS_BENCH_S").expect("RSOMICS_BENCH_S"));
    let g2m = PathBuf::from(std::env::var("RSOMICS_BENCH_G2M").expect("RSOMICS_BENCH_G2M"));
    let m = read_10x(&dir).unwrap();
    let s_genes = read_gene_list(&s).unwrap();
    let g2m_genes = read_gene_list(&g2m).unwrap();
    let params = CellCycleParams {
        n_bins: 25,
        seed: 0,
    };
    c.bench_function("cell_cycle", |b| {
        b.iter(|| black_box(score_cell_cycle(&m, &s_genes, &g2m_genes, &params).unwrap()))
    });
}

criterion_group!(benches, bench_cell_cycle);
criterion_main!(benches);
