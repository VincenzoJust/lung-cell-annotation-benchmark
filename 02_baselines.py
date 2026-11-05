"""PCA baseline embedding + CellTypist trained on the same reference."""
import celltypist
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from common import hvg, load

ref, qry = load()
genes = hvg(ref)

pca = PCA(50, random_state=0).fit(ref[:, genes].X.toarray())
emb = np.vstack([pca.transform(a[:, genes].X.toarray()) for a in (ref, qry)])
np.save("results/emb_pca.npy", emb.astype(np.float32))

# CellTypist densifies its input: all ~26k genes x 100k cells would need ~10 GB, so train on 5k HVGs
ct = hvg(ref, 5000)
# (the HVG subset no longer sums to 1e4 per cell, hence check_expression=False; annotate gets all genes)
model = celltypist.train(ref[:, ct], labels="coarse", use_SGD=True, mini_batch=True, n_jobs=-1, check_expression=False)
res = celltypist.annotate(qry, model=model)
pd.DataFrame({"pred": res.predicted_labels["predicted_labels"], "conf": res.probability_matrix.max(axis=1)}).to_csv("results/pred_celltypist.csv")
