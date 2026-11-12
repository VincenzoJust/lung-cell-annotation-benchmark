# Can AI models label lung cells in COVID-19 patients?

## The idea in plain words

Every cell in the body has a type (immune cell, airway cell, fibroblast...). Researchers measure which genes each cell uses and then label it, usually by comparing it to a large, well-annotated healthy "atlas".

This project asks: **when the cells come from sick lungs, which method labels them best?**

- **Atlas (known labels):** Human Lung Cell Atlas, healthy lungs (96,364 cells).
- **Patient cells (to label):** 116,313 cell nuclei from 7 control and 20 COVID-19 autopsy lungs (Melms et al. 2021).
- **Methods compared:** PCA, CellTypist, scVI/scANVI, and two "foundation models" (Geneformer, scGPT) that are pretrained on millions of cells and used as-is.
- **Score:** macro-F1 against the labels the original authors gave (1.0 = perfect, 0 = useless).

## Result

**Foundation models come close, but don't win.** Zero-shot Geneformer ties the best classic method (scANVI) on COVID cells (0.731 vs 0.734) without ever seeing them. All methods lose a similar amount going from control to COVID cells, so none is clearly more robust to disease. scGPT also invents 822 ionocytes (a rare cell type) that the authors never found.

| Method | Control | COVID-19 | Drop |
|---|---|---|---|
| scANVI (own classifier) | 0.812 | **0.734** | 0.078 |
| Geneformer + kNN | 0.807 | 0.731 | 0.077 |
| scANVI + kNN | 0.788 | 0.716 | 0.072 |
| scGPT + kNN | 0.797 | 0.711 | 0.086 |
| scVI + kNN | 0.780 | 0.681 | 0.099 |
| PCA + kNN | 0.506 | 0.495 | 0.010 |
| CellTypist | 0.463 | 0.438 | 0.025 |

(kNN = label a cell by its 15 nearest atlas cells.)

![label transfer](figures/label_transfer_f1.png)

Other findings, with plots in [figures/](figures/):
- **Integration:** scANVI mixes datasets best (scib-metrics total 0.659); the foundation models keep biology but don't remove batch effects.
- **COVID-specific states:** a collagen-rich fibroblast cluster (COL1A1/COL3A1/CTHRC1) and two macrophage clusters are more abundant in COVID. The models map them onto healthy types with high confidence, so low confidence did not flag them as new.
- **CFTR:** no difference between COVID and control in any cell type.
- **Differential expression:** thousands of genes differ, but the top hits (ribosomal, heat-shock) look like autopsy artifacts more than disease.

## Caveats

- One cohort of autopsy tissue; controls died of other causes.
- Atlas = whole cells, patients = nuclei. This technology gap is mixed up with the disease effect.
- Coarse labels only (15 lineages), one run, one seed, no confidence intervals.
- Foundation models are zero-shot only (no fine-tuning, 8 GB GPU).

## What's in this folder

| File | What it does |
|---|---|
| `01_get_data.py` | Downloads both datasets from CELLxGENE (~6.4 GB) |
| `02_baselines.py` | PCA and CellTypist |
| `03_scarches.py` | scVI / scANVI |
| `04_embed_fm.py` | Geneformer and scGPT |
| `analysis.ipynb` | Scores everything, makes all figures (outputs already included) |
| `label_map.csv` | Maps both datasets' cell-type names to 15 shared lineages |
| `common.py` | Small shared helpers |
| `figures/` | Plots used above |

Running the scripts creates `data/`, `models/`, `results/` and `logs/` (not in the repo).

## Reproduce

Needs Windows or Linux, Python 3.11 and a GPU (took ~6 h on an RTX 2060; Geneformer is ~3 h of that). Two environments are needed because the foundation-model package pins its own library versions. Install GPU torch first (see headers of the `requirements` files).

```bash
uv venv .venv --python 3.11
uv pip install --python .venv -r requirements.txt
uv venv .venv-fm --python 3.11
uv pip install --python .venv-fm -r requirements-fm.txt

.venv/Scripts/python 01_get_data.py 50     # optional smoke test on a tiny subset
.venv/Scripts/python 01_get_data.py
.venv/Scripts/python 02_baselines.py
.venv/Scripts/python 03_scarches.py
.venv-fm/Scripts/python 04_embed_fm.py scgpt
.venv-fm/Scripts/python 04_embed_fm.py geneformer
.venv/Scripts/jupyter nbconvert --to notebook --execute analysis.ipynb --inplace
```

On Linux/macOS use `.venv/bin/python` instead of `.venv/Scripts/python`.
