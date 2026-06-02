#!/usr/bin/env python3
"""scanpy oracle for rsomics-sc-cell-cycle.

Reads a 10x MTX directory (var_names='gene_ids' so identity matches the tool,
which keys on features.tsv column 1), runs sc.tl.score_genes_cell_cycle with
the given S and G2M gene lists (random_state=0 default => deterministic control
sampling), and writes a cell_id<TAB>S_score<TAB>G2M_score<TAB>phase TSV.

Usage: scanpy_cell_cycle_oracle.py <mtx_dir> <s_genes.txt> <g2m_genes.txt>
                                   <out.tsv> [n_bins seed]
"""
import sys

import scanpy as sc


def read_list(path):
    with open(path) as f:
        return [ln.strip() for ln in f if ln.strip()]


def main() -> None:
    mtx_dir, s_path, g2m_path, out_path = sys.argv[1:5]
    n_bins = int(sys.argv[5]) if len(sys.argv) > 5 else 25
    seed = int(sys.argv[6]) if len(sys.argv) > 6 else 0

    s_genes = read_list(s_path)
    g2m_genes = read_list(g2m_path)

    adata = sc.read_10x_mtx(mtx_dir, var_names="gene_ids")
    sc.tl.score_genes_cell_cycle(
        adata,
        s_genes=s_genes,
        g2m_genes=g2m_genes,
        n_bins=n_bins,
        random_state=seed,
    )

    with open(out_path, "w") as f:
        f.write("cell_id\tS_score\tG2M_score\tphase\n")
        for bc, s, g, p in zip(
            adata.obs_names,
            adata.obs["S_score"],
            adata.obs["G2M_score"],
            adata.obs["phase"],
        ):
            f.write(f"{bc}\t{float(s)!r}\t{float(g)!r}\t{p}\n")


if __name__ == "__main__":
    main()
