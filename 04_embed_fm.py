"""Zero-shot cell embeddings from a foundation model (via helical), ref rows then query rows.

.venv-fm/Scripts/python 04_embed_fm.py geneformer
.venv-fm/Scripts/python 04_embed_fm.py scgpt
"""
import sys

import anndata as ad
import numpy as np
import pandas as pd

name = sys.argv[1]
if name == "geneformer":
    from helical.models.geneformer import Geneformer, GeneformerConfig
    model = Geneformer(GeneformerConfig(model_name="gf-12L-104M-i4096", batch_size=16, device="cuda"))
else:
    from helical.models.scgpt import scGPT, scGPTConfig
    model = scGPT(scGPTConfig(batch_size=32, device="cuda"))

embs = []
for split in ("ref", "query"):
    a = ad.read_h5ad(f"data/{split}.h5ad")
    a.var["ensembl_id"] = a.var_names  # Geneformer: Ensembl ids; scGPT: symbols
    if name == "geneformer":  # keep 1:1 mappable genes (helical's duplicate-collapse path drops var['ensembl_id'])
        ids = pd.Series(a.var_names.str.upper()).map(model.tk.gene_mapping_dict)
        a = a[:, (ids.notna() & ~ids.duplicated()).values].copy()
    genes = "ensembl_id" if name == "geneformer" else "feature_name"
    embs.append(model.get_embeddings(model.process_data(a, gene_names=genes)))
    print(split, embs[-1].shape)

np.save(f"results/emb_{name}.npy", np.vstack(embs).astype(np.float32))
