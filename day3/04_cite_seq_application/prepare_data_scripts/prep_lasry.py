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
#     "seaborn>=0.13.2"
# ]
# ///


"""
Assemble Lasry 2022 CITE-seq RNA data from GEO GSE185381_RAW into a single AnnData.
"""

import argparse
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad


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
    # GSM5613768 = scRNA-seq only, skipped
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Assemble Lasry 2022 CITE-seq RNA AnnData from GEO GSE185381_RAW."
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
    return parser.parse_args()


def load_one_gsm(gsm: str, gse_dir: Path) -> ad.AnnData:
    """Load a single GSM as RNA-only AnnData with metadata merged."""
    cur_files = [x.name for x in gse_dir.iterdir() if gsm in x.name]

    def pick_one(pattern: str) -> str:
        matches = [f for f in cur_files if pattern in f]
        if len(matches) != 1:
            raise RuntimeError(f"Expected 1 match for pattern {pattern!r} in {gsm}, got {len(matches)}")
        return matches[0]

    matrix_file = pick_one("matrix.mtx.gz")
    features_file = pick_one("_features.tsv.gz")
    barcodes_file = pick_one("barcodes.tsv.gz")
    metadata_file = pick_one("metadata.csv.gz")

    # Read matrix and transpose to cells × genes
    ad_tmp = sc.read_mtx(gse_dir / matrix_file).T

    # Features
    features = pd.read_table(gse_dir / features_file, header=None)
    features.columns = ["gene_id", "symbol", "assay"]
    features.index = features["symbol"].to_numpy()
    ad_tmp.var = features

    # Barcodes
    barcodes = pd.read_table(gse_dir / barcodes_file, header=None)
    ad_tmp.obs_names = barcodes.iloc[:, 0].values

    # Keep only RNA (drop Antibody Capture)
    ad_rna = ad_tmp[:, ad_tmp.var["assay"] != "Antibody Capture"].copy()

    # Metadata
    meta = pd.read_csv(gse_dir / metadata_file, header=0)

    # Check that metadata barcodes match adata barcodes (ignoring 10x suffix)
    meta_barcodes = meta["cell"].str.split(":").apply(lambda x: x[1]).to_numpy()
    adata_barcodes = (
        pd.Series(ad_rna.obs_names)
        .str.split("-")
        .apply(lambda x: x[0])
        .to_numpy()
    )

    meta_barcodes_not_in_adata = (~np.isin(meta_barcodes, adata_barcodes)).sum()
    assert meta_barcodes_not_in_adata == 0, f"Barcode mismatch for {gsm}"

    # Align and merge metadata into obs
    mask = np.isin(adata_barcodes, meta_barcodes)
    ad_rna = ad_rna[mask, :].copy()
    assert ad_rna.n_obs == meta.shape[0], f"Cell count mismatch after masking for {gsm}"

    meta["barcode"] = meta["cell"].str.split(":").apply(lambda x: x[1])
    ad_rna.obs["barcode"] = (
        pd.Series(ad_rna.obs_names)
        .str.split("-")
        .apply(lambda x: x[0])
        .to_numpy()
    )

    ad_rna.obs = (
        ad_rna.obs.reset_index(names=["index"])
        .merge(meta, on="barcode", how="left")
        .set_index("index")
    )

    # annotate sample id
    ad_rna.obs["sample"] = gsm

    return ad_rna


def build_lasry_reference(gse_dir: Path) -> ad.AnnData:
    """Load all CITE-seq GSMs and concatenate into a single AnnData."""
    if not gse_dir.exists():
        raise FileNotFoundError(f"GSE directory not found: {gse_dir}")

    rna_adatas: Dict[str, ad.AnnData] = {}

    for gsm in CITE_GSMS:
        print(f"Processing {gsm} ...")
        ad_rna = load_one_gsm(gsm, gse_dir)
        ad_rna.var_names_make_unique()
        rna_adatas[gsm] = ad_rna

    # concatenate along cells, keeping sample label
    adata = ad.concat(
        rna_adatas,
        label="sample_from_concat",
        join="outer",
        merge="unique",
        index_unique="-",
    )

    # Keep a clean "sample" column (we already set obs['sample'] above)
    # You can drop sample_from_concat if redundant:
    # del adata.obs["sample_from_concat"]

    return adata


def main() -> None:
    args = parse_args()
    gse_dir: Path = args.gse_dir
    output_path: Path = args.output

    print(f"GSE directory: {gse_dir}")
    print(f"Output file:   {output_path}")

    adata = build_lasry_reference(gse_dir)

    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Writing AnnData to {output_path}")
    import anndata
    anndata.settings.allow_write_nullable_strings = True
    adata.write_h5ad(output_path)

    print("Done.")


if __name__ == "__main__":
    main()

