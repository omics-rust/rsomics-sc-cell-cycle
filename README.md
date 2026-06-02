# rsomics-sc-cell-cycle

Per-cell cell-cycle phase scoring for single-cell expression. Given an S-phase
and a G2M-phase gene list, it computes an `S_score` and a `G2M_score` per cell
(each a `score_genes`-style set mean minus a bin-matched control mean) and
assigns each cell a phase — `S`, `G2M`, or `G1`. This is scanpy's
`sc.tl.score_genes_cell_cycle`.

```
rsomics-sc-cell-cycle <10x-mtx-dir> --s-genes s.txt --g2m-genes g2m.txt \
    -o phase.tsv [--n-bins 25 --seed 0]
```

Input is a 10x MatrixMarket directory (`matrix.mtx[.gz]`, `features.tsv[.gz]`,
`barcodes.tsv[.gz]`) plus two gene-list files (one gene id — `features.tsv`
column 1 — per line). Output is a `cell_id<TAB>S_score<TAB>G2M_score<TAB>phase`
TSV, one line per cell.

The two scores are computed with `ctrl_size = min(len(s_genes),
len(g2m_genes))`, exactly as scanpy fixes it. The phase call is then: `G1` if
both scores are negative, else `G2M` if the G2M score exceeds the S score, else
`S`. With `--seed 0` the control draw — and therefore every score and phase — is
value-exact against `scanpy.tl.score_genes_cell_cycle(adata, random_state=0)`.

## Origin

Independent Rust reimplementation of scanpy's `sc.tl.score_genes_cell_cycle`,
based on the scanpy source (`scanpy/tools/_score_genes.py`, BSD-3-Clause —
readable, cited) and the Seurat cell-cycle scoring method:

- Satija et al., "Spatial reconstruction of single-cell gene expression data",
  Nature Biotechnology 33, 495–502 (2015), DOI 10.1038/nbt.3192.

The scoring core (bin-matched control sampling reproducing numpy's `RandomState`
MT19937 + Fisher-Yates `choice(replace=False)` bit-for-bit) is shared with
`rsomics-sc-score-genes`; this crate adds the deterministic two-score phase call.

License: MIT OR Apache-2.0. Upstream credit: scanpy
(https://github.com/scverse/scanpy, BSD-3-Clause).
