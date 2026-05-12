import marimo

__generated_with = "0.23.5"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Single-cell multi-omics analysis of the immune response in COVID-19
    Wellcome HCA Strategic Science Support
    Analysis of human blood immune cells provides insights into the coordinated response to viral infections such as severe acute respiratory syndrome coronavirus 2, which causes coronavirus disease 2019 (COVID-19). We performed single-cell transcriptome, surface proteome and T and B lymphocyte antigen receptor analyses of over 780,000 peripheral blood mononuclear cells from a cross-sectional cohort of 130 patients with varying severities of COVID-19. We identified expansion of nonclassical monocytes expressing complement transcripts (CD16+C1QA/B/C+) that sequester platelets and were predicted to replenish the alveolar macrophage pool in COVID-19. Early, uncommitted CD34+ hematopoietic stem/progenitor cells were primed toward megakaryopoiesis, accompanied by expanded megakaryocyte-committed progenitors and increased platelet activation. Clonally expanded CD8+ T cells and an increased ratio of CD8+ effector T cells to effector memory T cells characterized severe disease, while circulating follicular helper T cells accompanied mild disease. We observed a relative loss of IgA2 in symptomatic disease despite an overall expansion of plasmablasts and plasma cells. Our study highlights the coordinated immune response that contributes to COVID-19 pathogenesis and reveals discrete cellular components that can be targeted for therapy.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    https://www.nature.com/articles/s41591-021-01329-2
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    https://cellxgene.cziscience.com/collections/ddfad306-714d-4cc0-9985-d9072820c530
    """)
    return


@app.cell
def _():
    import scanpy as sc
    import polars as pl
    import os

    return os, sc


@app.cell
def _(os):
    data_dir = "/home/rstudio/project/XX_exam_data/data"
    input_dir = os.path.join(data_dir, "haniffa_stevenson_covid19.h5ad")
    return data_dir, input_dir


@app.cell
def _(input_dir, sc):
    adata = sc.read_h5ad(input_dir)
    return (adata,)


@app.cell
def _(adata):
    adata
    return


@app.cell
def _(adata):
    import numpy as np
    import scipy.sparse as sp

    X = adata.X
    total_counts = adata.obs["total_counts"].to_numpy()
    target_sum = 1e4

    if sp.issparse(X):
        X_norm = X.copy()
        X_norm.data = np.expm1(X_norm.data)
        raw_approx = X_norm.multiply(total_counts[:, None] / target_sum)
        raw_approx.data = np.rint(raw_approx.data)
    else:
        X_norm = np.expm1(X)
        raw_approx = np.rint(X_norm * (total_counts[:, None] / target_sum)).astype(int)
    return np, raw_approx, sp


@app.cell
def _(adata, raw_approx):
    adata.layers["counts"] = raw_approx.tocsr()
    return


@app.cell
def _(adata):
    adata.uns["antibody_raw.X"].shape
    return


@app.cell
def _(np, sp):
    import pandas as pd
    from anndata import AnnData
    from geosketch import gs

    def geosketch_rna_adt_from_uns(adata, n_cells: int=10000, rep: str='X_pca', seed: int=1, adt_counts_key: str='antibody_raw.X', adt_features_key: str='antibody_features', sketch_key: str='is_geosketch'):
        """
        Geosketch cells from an AnnData object using an embedding in adata.obsm,
        while manually subsetting ADT counts stored in adata.uns.

        Returns
        -------
        rna_sketch
            AnnData object with RNA data.
        adt_sketch
            AnnData object with ADT counts.
        sketch_idx
            Integer indices of selected cells.
        """
        if rep not in adata.obsm:
            raise KeyError(f'{rep!r} not found in adata.obsm')
        if adt_counts_key not in adata.uns:
            raise KeyError(f'{adt_counts_key!r} not found in adata.uns')
        if adt_features_key not in adata.uns:
            raise KeyError(f'{adt_features_key!r} not found in adata.uns')
        if n_cells > adata.n_obs:
            raise ValueError(f'n_cells={n_cells} is larger than adata.n_obs={adata.n_obs}')
        np.random.seed(seed)
        X_rep = adata.obsm[rep]
        sketch_idx = gs(X_rep, n_cells, replace=False)
        sketch_idx = np.asarray(sketch_idx)
        sketch_idx = np.sort(sketch_idx)
        adata.obs[sketch_key] = False
        adata.obs.loc[adata.obs_names[sketch_idx], sketch_key] = True
        rna_sketch = adata[sketch_idx].copy()
        for key in [adt_counts_key, 'antibody_X', adt_features_key]:
            if key in rna_sketch.uns:
                del rna_sketch.uns[key]
        adt_counts = adata.uns[adt_counts_key]
        if sp.issparse(adt_counts):
            adt_counts_sketch = adt_counts[sketch_idx, :].copy()
        else:
            adt_counts_sketch = np.asarray(adt_counts)[sketch_idx, :].copy()
        adt_features = adata.uns[adt_features_key]
        if isinstance(adt_features, pd.DataFrame):
            adt_var = adt_features.copy()
            if adt_var.index.size != adt_counts_sketch.shape[1]:
                raise ValueError('antibody_features has a different number of rows than ADT features.')  # 1. Select representative cells using RNA PCA or another embedding
        else:
            adt_names = pd.Index(np.asarray(adt_features).astype(str), name='feature_name')
            adt_var = pd.DataFrame(index=adt_names)
        adt_var.index = adt_var.index.astype(str)  # Make ordering stable/reproducible in the output object
        adt_var_names = pd.Index(adt_var.index)
        if not adt_var_names.is_unique:
            adt_var.index = pd.Index([f'{name}_{i}' for i, name in enumerate(adt_var.index)], name=adt_var.index.name)
        adt_obs = rna_sketch.obs.copy()  # 2. Mark selected cells in the full object
        adt_sketch = AnnData(X=adt_counts_sketch, obs=adt_obs, var=adt_var)
        adt_sketch.var['feature_type'] = 'Antibody Capture'
        return (rna_sketch, adt_sketch, sketch_idx)  # 3. Subset RNA normally  # Optional: remove the bulky ADT objects from RNA .uns  # because we will save ADT as its own object.  # 4. Manually subset ADT counts stored in .uns  # 5. Build ADT feature metadata  # Ensure ADT var names are unique  # 6. Use same obs as RNA sketch

    return (geosketch_rna_adt_from_uns,)


@app.cell
def _(adata, geosketch_rna_adt_from_uns):
    rna_sketch, adt_sketch, sketch_idx = geosketch_rna_adt_from_uns(
        adata,
        n_cells=30_000,
        rep="X_pca",
        seed=1,
    )
    return adt_sketch, rna_sketch


@app.cell
def _(rna_sketch):
    rna_sketch
    return


@app.cell
def _(adt_sketch):
    adt_sketch.var
    return


@app.cell
def _():
    import pandas as pd
    import anndata as ad

    ad.settings.allow_write_nullable_strings = True

    return


@app.cell
def _(adt_sketch, data_dir, os, rna_sketch):
    rna_sketch.write_h5ad(
        os.path.join(data_dir, "exam_rna_geosketch_30k.h5ad"),
        compression="gzip",
    )

    adt_sketch.write_h5ad(
        os.path.join(data_dir, "exam_adt_geosketch_30k.h5ad"),
        compression="gzip",
    )
    return


if __name__ == "__main__":
    app.run()
