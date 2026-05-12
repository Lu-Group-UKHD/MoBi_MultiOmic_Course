# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "igraph>=1.0.0",
#     "ipykernel>=7.2.0",
#     "ipywidgets>=8.1.8",
#     "jupyterlab>=4.5.7",
#     "matplotlib>=3.10.9",
#     "nbconvert>=7.17.1",
#     "numpy>=2.4.4",
#     "pandas>=3.0.2",
#     "scanpy>=1.12.1",
#     "seaborn>=0.13.2",
#     "mudata",
#     "muon",
#     "geosketch"
# ]
# ///

#!/usr/bin/env python3
"""
Assemble Lasry 2022 CITE-seq data from GEO GSE185381_RAW.

Per GSM:
- read raw 10x-like matrix
- assign features + barcodes
- keep RNA rows only
- subset to metadata-defined cells
- attach processed ADT if available
- return MuData({"rna": ..., "adt": ...}) or RNA-only MuData

Finally:
- concatenate RNA across all GSMs
- save merged RNA .h5ad
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from geosketch import gs



# PROJECT_ROOT = Path("/home/rstudio/project/aml-lsc-analysis")
# GSE_DIR = PROJECT_ROOT / "data" / "GSE185381_RAW"
# OUT_DIR = PROJECT_ROOT / "data_clean" / "assembled_artifacts" / "gse185381_lasry_citeseq"
# OUT_DIR.mkdir(parents=True, exist_ok=True)

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Assemble Lasry 2022 CITE-seq RNA+ADT AnnData from GEO GSE185381_RAW."
    )
    parser.add_argument(
        "--gse-dir",
        type=Path,
        default=Path.home() / "data/aml-lsc/raw/lasry_2022/GSE185381_RAW",
        help="Directory containing GSE185381_RAW files (matrix.mtx.gz, features.tsv.gz, barcodes, metadata). "
             "Default: ~/data/aml-lsc/raw/lasry_2022/GSE185381_RAW",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path.home() / "data/aml-lsc/processed/lasry_2022/lasry_rna_assembled.h5ad",
        help="Output .h5ad file path. Default: ~/data/aml-lsc/processed/lasry_2022/lasry_rna_assembled.h5ad",
    )
    parser.add_argument(
        "--n_sketch",
        type=int,
        default=30000,
        help="Number of cells to include in the sketch. Default: 30000",
    )
    return parser.parse_args()

args = parse_args()
GSE_DIR =  args.gse_dir
OUT_DIR = args.output.parent
N_SKETCH = args.n_sketch

CITE_GSMS = [
    "GSM5613744", "GSM5613745", "GSM5613746",
    "GSM5613747", "GSM5613748", "GSM5613749", "GSM5613750",
    "GSM5613751", "GSM5613752",
    "GSM5613753", "GSM5613754", "GSM5613755",
    "GSM5613756",
    "GSM5613757", "GSM5613758", "GSM5613759",
    "GSM5613760", "GSM5613761",
    "GSM5613762", "GSM5613763",
    "GSM5613764", "GSM5613765",
    "GSM5613766", "GSM5613767",
    # GSM5613768 is scRNA-only, skipped
    "GSM5613769", "GSM5613770", "GSM5613771",
    "GSM5613772", "GSM5613773",
    "GSM5613774", "GSM5613775", "GSM5613776",
    "GSM5613777", "GSM5613778",
    "GSM5613779", "GSM5613780",
    "GSM5613781",
    "GSM5613782",
    "GSM5613783", "GSM5613784", "GSM5613785", "GSM5613786",
    "GSM5613787",
    "GSM5613788", "GSM5613789", "GSM5613790",
]


def log(*args):
    print(*args, flush=True)


def list_gsm_files(gsm: str, gse_dir: Path) -> list[Path]:
    return sorted(
        p for p in gse_dir.iterdir()
        if p.is_file() and p.name.startswith(gsm)
    )


def pick_one(files: list[Path], pattern: str, required: bool = True) -> Path | None:
    hits = [f for f in files if pattern in f.name]
    if len(hits) == 1:
        return hits[0]
    if len(hits) == 0:
        if required:
            raise RuntimeError(f"Missing file matching {pattern!r}")
        return None
    raise RuntimeError(f"Expected 1 file matching {pattern!r}, found {[x.name for x in hits]}")


def parse_sample_key(gsm: str, matrix_file: Path) -> str:
    m = re.match(rf"^{gsm}_(.+)_matrix\.mtx\.gz$", matrix_file.name)
    return m.group(1) if m else gsm


def clean_obs_names(index_like) -> pd.Index:
    idx = pd.Index(pd.Series(index_like, dtype="string").astype(str).to_numpy(), dtype="object")
    idx.name = None
    return idx


def extract_base_barcodes(barcodes) -> pd.Index:
    """
    Turn 10x-style cell barcodes like AAAC...-1 into AAAC...
    """
    s = pd.Series(barcodes, dtype="string").astype(str)
    out = s.str.split("-").str[0]
    idx = pd.Index(out.to_numpy(), dtype="object")
    idx.name = None
    return idx


def read_features(features_file: Path) -> pd.DataFrame:
    features = pd.read_table(features_file, header=None)

    if features.shape[1] == 3:
        features.columns = ["gene_id", "symbol", "assay"]
    elif features.shape[1] == 2:
        features.columns = ["gene_id", "symbol"]
        features["assay"] = "Gene Expression"
    else:
        raise RuntimeError(f"Unexpected features shape for {features_file.name}: {features.shape}")

    features["gene_id"] = features["gene_id"].astype(str)
    features["symbol"] = features["symbol"].astype(str)
    features["assay"] = features["assay"].astype(str)

    features.index = pd.Index(features["symbol"].to_numpy(), dtype="object")
    features.index.name = None
    return features


def read_metadata(metadata_file: Path) -> pd.DataFrame:
    meta = pd.read_csv(metadata_file)

    if "cell" not in meta.columns:
        raise RuntimeError(f"'cell' column not found in {metadata_file.name}")

    # robustly grab barcode after the colon
    meta["barcode"] = meta["cell"].astype(str).str.split(":").str[-1]
    meta["barcode"] = meta["barcode"].astype(str)

    # keep one row per barcode
    meta = meta.drop_duplicates("barcode").copy()
    return meta


def load_adt(adt_file: Path, sample_key: str, rna_base_barcodes: pd.Index, obs_index: pd.Index) -> ad.AnnData:
    """
    ADT_processed.csv.gz is features x cells, dense.
    We transpose to cells x features and align to RNA cell order.
    """
    adt = pd.read_csv(adt_file, index_col=0)

    adt.columns = (
        pd.Series(adt.columns, dtype="string")
        .astype(str)
        .str.replace(f"^{re.escape(sample_key.replace('-', '.'))}\\.", "", regex=True)
        .str.replace(".", "-", regex=False)
        .to_numpy()
    )

    adt = adt.T
    adt.index = clean_obs_names(adt.index)

    adt_base = extract_base_barcodes(adt.index)
    adt["__barcode__"] = adt_base.to_numpy()
    adt = adt.drop_duplicates("__barcode__").set_index("__barcode__")

    adt = adt.reindex(rna_base_barcodes)
    X = sp.csr_matrix(adt.fillna(0).to_numpy())

    adata_adt = ad.AnnData(X=X)
    adata_adt.obs_names = obs_index.copy()
    adata_adt.var_names = clean_obs_names(adt.columns)
    adata_adt.var["feature_type"] = "Antibody Capture"

    return adata_adt


def load_one_gsm(gsm: str, gse_dir: Path = GSE_DIR) -> md.MuData:
    files = list_gsm_files(gsm, gse_dir)

    matrix_file = pick_one(files, "matrix.mtx.gz")
    features_file = pick_one(files, "_features.tsv.gz")
    barcodes_file = pick_one(files, "barcodes.tsv.gz")
    metadata_file = pick_one(files, "metadata.csv.gz")
    adt_file = pick_one(files, "ADT_processed.csv.gz", required=False)

    sample_key = parse_sample_key(gsm, matrix_file)
    log(f"Loading {gsm} ({sample_key})")

    # Read matrix and transpose to cells x features
    ad_tmp = sc.read_mtx(matrix_file).T

    # Features
    features = read_features(features_file)
    ad_tmp.var = features
    ad_tmp.var_names = clean_obs_names(features.index)
    ad_tmp.var_names_make_unique()

    # Barcodes
    barcodes = pd.read_table(barcodes_file, header=None).iloc[:, 0].astype(str)
    ad_tmp.obs_names = clean_obs_names(barcodes)

    # RNA only
    ad_rna = ad_tmp[:, ad_tmp.var["assay"].ne("Antibody Capture").to_numpy()].copy()

    # Metadata and filtering to called cells
    meta = read_metadata(metadata_file)
    meta_barcodes = pd.Index(meta["barcode"].astype(str).to_numpy(), dtype="object")
    meta_barcodes.name = None

    rna_base_barcodes = extract_base_barcodes(ad_rna.obs_names)
    keep = rna_base_barcodes.isin(meta_barcodes)

    ad_rna = ad_rna[keep, :].copy()
    rna_base_barcodes = extract_base_barcodes(ad_rna.obs_names)

    # Align metadata to RNA cell order
    meta_aligned = (
        meta.set_index("barcode")
        .reindex(rna_base_barcodes)
        .reset_index()
    )

    # Create unique obs names by appending sample key
    final_obs_names = clean_obs_names([f"{bc}-{sample_key}" for bc in rna_base_barcodes])
    ad_rna.obs_names = final_obs_names

    # Build obs cleanly
    obs = pd.DataFrame(index=final_obs_names)
    obs.index.name = None
    obs["barcode"] = rna_base_barcodes.to_numpy()
    obs["gsm_id"] = gsm
    obs["sample_key"] = sample_key
    obs["sample"] = gsm
    obs["dataset"] = "GSE185381"
    obs["source_dataset"] = "GSE185381"
    obs["modality"] = "CITE-seq"

    for col in meta_aligned.columns:
        if col == "barcode":
            continue
        obs[col] = meta_aligned[col].to_numpy()

    if "samples" in obs.columns:
        obs["patient_id"] = obs["samples"]
        obs["aml_id"] = obs["samples"]

    ad_rna.obs = obs

    mods = {"rna": ad_rna}

    if adt_file is not None:
        try:
            ad_adt = load_adt(
                adt_file=adt_file,
                sample_key=sample_key,
                rna_base_barcodes=rna_base_barcodes,
                obs_index=final_obs_names,
            )
            ad_adt.obs = obs.copy()
            mods["adt"] = ad_adt
        except Exception as e:
            log(f"  WARNING: failed to attach ADT for {gsm}: {e}")

    mdata = md.MuData(mods)
    return mdata


def build_all(gse_dir: Path = GSE_DIR) -> tuple[dict[str, md.MuData], ad.AnnData]:
    mdata_dict: dict[str, md.MuData] = {}

    for gsm in CITE_GSMS:
        try:
            mdata = load_one_gsm(gsm, gse_dir)
            print(mdata)
            print(mdata["adt"].X.data)
            mdata_dict[gsm] = mdata
            log(
                f"  -> RNA: {mdata['rna'].n_obs} cells x {mdata['rna'].n_vars} genes"
                + (f" | ADT: {mdata['adt'].n_obs} x {mdata['adt'].n_vars}" if "adt" in mdata.mod else "")
            )
        except Exception as e:
            log(f"ERROR loading {gsm}: {e}")

    if not mdata_dict:
        raise RuntimeError("No samples loaded")

    rna_list = [m["rna"] for m in mdata_dict.values()]
    merged_rna = ad.concat(
        rna_list,
        axis=0,
        join="outer",
        merge="unique",
        fill_value=0,
    )
    
    adt_list = [m["adt"] for m in mdata_dict.values()]
    merged_adt = ad.concat(
        adt_list,
        axis=0,
        join="outer",
        merge="unique",
        fill_value=0,
    )

    return mdata_dict, merged_rna, merged_adt


def main():
    log("=== GSE185381 Lasry assembly ===")
    mdata_dict, merged_rna, merged_adt = build_all(GSE_DIR)

    log(f"Merged RNA: {merged_rna.n_obs} cells x {merged_rna.n_vars} genes")
    log(f"Merged ADT: {merged_adt.n_obs} cells x {merged_adt.n_vars} features")

    detected = np.asarray((merged_rna.X > 0).sum(axis=1)).ravel()
    log("Cells with 0 genes:", int((detected == 0).sum()))
    log("Cells with <200 genes:", int((detected < 200).sum()))

    # out_file = OUT_DIR / "gse185381_lasry_rna.h5ad"
    # merged_rna.write_h5ad(out_file)
    
    import anndata
    anndata.settings.allow_write_nullable_strings = True
    # mdata.write_h5mu(args.output)
    
    merged_rna.write_h5ad(args.output)
    merged_adt = merged_adt[merged_rna.obs_names, :].copy()
    
    merged_adt.write_h5ad(str(args.output).replace(".h5ad", "_adt.h5ad"))
    
    log("Saved:", args.output)
    
    
    bdata = sc.AnnData(
        X = merged_rna.X.copy(),
        obs = merged_rna.obs.copy(),
        var=merged_rna.var.copy(),
    )
    sc.pp.normalize_total(bdata, target_sum=1e4)
    sc.pp.log1p(bdata)
    sc.pp.highly_variable_genes(bdata, n_top_genes=3000, batch_key="samples")

    sc.pp.pca(bdata, n_comps=50, use_highly_variable=True)
    
    mat = bdata.obsm["X_pca"]
    sketch_index = gs(mat, N=N_SKETCH, replace=False)
    
    rna_sketch = merged_rna[sketch_index, :].copy()
    adt_sketch = merged_adt[sketch_index, :].copy()
    
    rna_sketch.write_h5ad(str(args.output).replace(".h5ad", f"_rna_sketch_{N_SKETCH}.h5ad"))
    adt_sketch.write_h5ad(str(args.output).replace(".h5ad", f"_adt_sketch_{N_SKETCH}.h5ad"))
    
    log("Saved sketch:", str(args.output).replace(".h5ad", f"_rna_sketch_{N_SKETCH}.h5ad"))
    log("Saved sketch:", str(args.output).replace(".h5ad", f"_adt_sketch_{N_SKETCH}.h5ad"))

if __name__ == "__main__":
    main()
