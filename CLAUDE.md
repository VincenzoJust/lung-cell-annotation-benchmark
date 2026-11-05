# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A benchmark of cell-type annotation methods on diseased lung: the HLCA core (healthy reference, CELLxGENE) against Melms et al. 2021 COVID-19 lung (query: 7 control + 20 COVID donors, snRNA-seq). It compares PCA, CellTypist, scVI/scANVI (scArches surgery), Geneformer and scGPT (zero-shot), then covers disease-specific states, CFTR/ionocytes, and pseudobulk DE. Project scope/status lives in `BELLION/BERU/4_Templates/CV/Project-Portfolio.md`, not here.

## Environments (Windows, RTX 2060 8 GB, 16 GB RAM)

System Python is 3.14, which is too new. Two uv venvs on Python 3.11:

```bash
uv venv .venv --python 3.11
uv pip install --python .venv torch --index-url https://download.pytorch.org/whl/cu124
uv pip install --python .venv -r requirements.txt
uv venv .venv-fm --python 3.11
uv pip install --python .venv-fm -r requirements-fm.txt   # then reinstall GPU torch, see file header
```

`.venv` runs `01`–`03` and `analysis.ipynb`. `.venv-fm` runs `04` only, because helical pins its own torch and transformers.

## Pipeline

```bash
.venv/Scripts/python 01_get_data.py 50     # smoke test: 50 cells per type in ref AND query
.venv/Scripts/python 01_get_data.py        # full: <=3000 ref cells per cell type, all ~116k query cells
.venv/Scripts/python 02_baselines.py       # emb_pca + pred_celltypist
.venv/Scripts/python 03_scarches.py        # emb_scvi, emb_scanvi, pred_scanvi; models/ saved
.venv-fm/Scripts/python 04_embed_fm.py geneformer   # slow: ~25 cells/s
.venv-fm/Scripts/python 04_embed_fm.py scgpt        # ~130 cells/s
.venv/Scripts/jupyter nbconvert --to notebook --execute analysis.ipynb --inplace
```

Always run the smoke test end to end before a full run. On full data, run `03` and `04` in the background.

## Data contract between steps

- `01` writes `data/ref.h5ad` and `data/query.h5ad` with **raw counts in X**, Ensembl IDs as `var_names`, symbols in `var['feature_name']`, both restricted to the shared genes, and a `coarse` obs column. Ref cells with `coarse == "unmapped"` are dropped. Query cells keep it and are excluded only from scoring.
- Every `results/emb_<method>.npy` holds **ref rows first, then query rows**, in file order. The notebook asserts `len == n_ref + n_query` and fits one shared kNN (k=15) on the ref rows. A new method only has to write an embedding in this layout, or a `results/pred_<method>.csv` (index = query obs_names, columns `pred`, `conf`) and be added to the native-prediction list in the notebook.
- `common.load()` adds log1p-normalized X plus `layers['counts']`, and sets the batch column (`dataset` = HLCA study; query cells get their `donor_id`).

## Labels

`label_map.csv` maps **two different vocabularies** onto about 15 coarse lineages:
- reference: HLCA CL `cell_type`
- query: the authors' `cell_type_intermediate`. Never use the query's CL `cell_type`: CELLxGENE files Melms' author AT1 cells under "pulmonary alveolar type 2 cell".

Query pericytes and mesothelial cells sit under the authors' "Fibroblasts", so the shared class is `stromal`. `ionocyte` exists only in the reference. Quote labels that contain commas.

## Known gotchas

- CELLxGENE Census (`tiledbsoma`) has no Windows wheels, so `01` downloads h5ads from the CELLxGENE curation API with resumable `curl -C -` and reads the 5.9 GB HLCA in backed mode.
- seurat_v3 HVG with `batch_key` crashes on tiny per-study batches, so `common.hvg` runs it without a batch key.
- helical's Geneformer duplicate-collapse path loses `var['ensembl_id']`. `04` keeps only genes with a 1:1 mapping in `model.tk.gene_mapping_dict` to avoid that path.
- CellTypist densifies its input (all genes × 100k cells ≈ 10 GB), so `02` trains on 5k HVGs with `check_expression=False`, because the subset no longer sums to 1e4. `annotate` gets the full query.
- `pynndescent<0.7` is pinned: 0.7 returns `csr_array`, and scanpy 1.11 `neighbors` then fails on `getnnz`.
- `analysis.ipynb` writes nothing when a cell fails under `nbconvert --inplace`. Read `logs/analysis.log` for the error. The joblib `KeyError: ...joblib_memmapping_folder...` lines are harmless Windows cleanup noise.
- scVI is set up without `labels_key`. Labels are added in `SCANVI.from_scvi_model`, otherwise query surgery fails on the unseen "Unknown" category.

## Conventions

English. Plots use the Gold Tangerine palette from `BERU/4_Templates/Color-Templates/Color-Schemes.md` (already defined at the top of `analysis.ipynb`). `data/`, `models/`, `*.npy` and `*.h5ad` are gitignored.
