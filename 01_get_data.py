"""Download HLCA core + Melms 2021 COVID lung from CELLxGENE, subsample the reference, write raw-count h5ads.

python 01_get_data.py        # full: up to 3000 ref cells per cell type, all query cells
python 01_get_data.py 50     # smoke test: 50 per type in ref and query
"""
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import requests

API = "https://api.cellxgene.cziscience.com/curation/v1"
REF = ("6f6d381a-7701-4781-935c-db10d30de293", "066943a2-fdac-4b29-b348-40cede398e4e")  # HLCA core
QRY = ("e4c9ed14-e560-4900-a3bf-b0f8d2ce6a10", "d8da613f-e681-4c69-b463-e94f5e66847f")  # Melms 2021
CAP = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
SMALL = len(sys.argv) > 1
D = Path("data")
D.mkdir(exist_ok=True)


def download(coll, ds):
    out = D / f"{ds}.h5ad"
    if not out.exists():
        meta = requests.get(f"{API}/collections/{coll}/datasets/{ds}").json()
        url = next(a["url"] for a in meta["assets"] if a["filetype"] == "H5AD")
        print("downloading", meta["title"], url)
        # curl resumes (-C -) after dropped connections
        subprocess.run(["curl", "-L", "-sS", "-C", "-", "--retry", "20", "--retry-all-errors", "-o", f"{out}.part", url], check=True)
        Path(f"{out}.part").rename(out)
    return out


def stratified(obs, cap, seed=0):
    """Row positions, at most `cap` cells per cell_type (rare types kept whole)."""
    rng = np.random.default_rng(seed)
    keep = [rng.choice(g, min(cap, len(g)), replace=False) for g in obs.groupby("cell_type", observed=True).indices.values()]
    return np.sort(np.concatenate(keep))


def load(path, idx=None):
    """Raw counts in X, Ensembl ids as var_names, gene symbols in var['feature_name']."""
    a = ad.read_h5ad(path, backed="r")
    idx = np.arange(a.n_obs) if idx is None else idx
    raw = a.raw if a.raw is not None else a
    out = ad.AnnData(X=raw.X[idx], obs=a.obs.iloc[idx].copy(), var=raw.var[["feature_name"]].copy())
    out.var["feature_name"] = out.var["feature_name"].astype(str)
    a.file.close()
    return out


ref_path, qry_path = download(*REF), download(*QRY)

obs = ad.read_h5ad(ref_path, backed="r").obs
ref = load(ref_path, stratified(obs, CAP))
qry_obs = ad.read_h5ad(qry_path, backed="r").obs
qry = load(qry_path, stratified(qry_obs, CAP) if SMALL else None)

genes = ref.var_names.intersection(qry.var_names)
ref, qry = ref[:, genes].copy(), qry[:, genes].copy()

# Coarse labels shared by both datasets. Query truth comes from the authors' cell_type_intermediate:
# CELLxGENE's ontology mapping for Melms files author AT1 cells under "pulmonary alveolar type 2 cell".
lm = pd.read_csv("label_map.csv", index_col=0)["coarse"]
for a, col in ((ref, "cell_type"), (qry, "cell_type_intermediate")):
    a.obs["coarse"] = a.obs[col].astype(str).map(lm).fillna("unmapped")
    print(a.n_obs, "cells,", a.n_vars, "genes; unmapped:", sorted(a.obs.loc[a.obs.coarse == "unmapped", col].unique()))

ref = ref[ref.obs.coarse != "unmapped"].copy()  # reference needs a label
ref.write_h5ad(D / "ref.h5ad", compression="gzip")
qry.write_h5ad(D / "query.h5ad", compression="gzip")
