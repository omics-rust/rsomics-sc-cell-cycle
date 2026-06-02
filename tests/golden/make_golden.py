#!/usr/bin/env python3
"""Generate the committed golden: a small synthetic 10x dir, S and G2M gene-list
files, and the scanpy score_genes_cell_cycle reference TSV. Run once with the
scanpy venv; the outputs are committed so CI validates without scanpy installed.
"""
import gzip
import os
import sys

import numpy as np
import scanpy as sc

HERE = os.path.dirname(os.path.abspath(__file__))
TENX = os.path.join(HERE, "tenx")


def write_tenx(counts):
    os.makedirs(TENX, exist_ok=True)
    n_genes, n_cells = counts.shape
    rows, cols = np.nonzero(counts)
    with gzip.open(os.path.join(TENX, "matrix.mtx.gz"), "wt") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write(f"{n_genes} {n_cells} {rows.size}\n")
        for r, c in zip(rows, cols):
            f.write(f"{r + 1} {c + 1} {int(counts[r, c])}\n")
    with gzip.open(os.path.join(TENX, "features.tsv.gz"), "wt") as f:
        for g in range(n_genes):
            f.write(f"ENSG{g:08d}\tGene{g}\tGene Expression\n")
    with gzip.open(os.path.join(TENX, "barcodes.tsv.gz"), "wt") as f:
        for c in range(n_cells):
            f.write(f"CELL{c:08d}-1\n")


def main():
    rng = np.random.default_rng(11)
    n_genes, n_cells = 200, 80
    base = rng.uniform(0.1, 6.0, size=n_genes)
    counts = rng.poisson(base[:, None], size=(n_genes, n_cells))
    write_tenx(counts)

    gene_ids = [f"ENSG{g:08d}" for g in range(n_genes)]
    s_genes = gene_ids[10:30]
    g2m_genes = gene_ids[100:130]
    with open(os.path.join(HERE, "s_genes.txt"), "w") as f:
        f.write("\n".join(s_genes) + "\n")
    with open(os.path.join(HERE, "g2m_genes.txt"), "w") as f:
        f.write("\n".join(g2m_genes) + "\n")

    adata = sc.read_10x_mtx(TENX, var_names="gene_ids")
    sc.tl.score_genes_cell_cycle(
        adata, s_genes=s_genes, g2m_genes=g2m_genes, random_state=0
    )
    with open(os.path.join(HERE, "cell_cycle.tsv"), "w") as f:
        f.write("cell_id\tS_score\tG2M_score\tphase\n")
        for bc, s, g, p in zip(
            adata.obs_names,
            adata.obs["S_score"],
            adata.obs["G2M_score"],
            adata.obs["phase"],
        ):
            f.write(f"{bc}\t{float(s)!r}\t{float(g)!r}\t{p}\n")
    print(
        f"wrote golden: {n_genes} genes x {n_cells} cells, "
        f"S={len(s_genes)} G2M={len(g2m_genes)}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
