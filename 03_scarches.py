"""scVI and scANVI on the reference, scArches surgery on the query (unlabeled)."""
import numpy as np
import pandas as pd
import scvi

from common import BATCH, hvg, load

scvi.settings.seed = 0
ref, qry = load()
genes = hvg(ref)
ref, qry = ref[:, genes].copy(), qry[:, genes].copy()
qry.obs["coarse_truth"] = qry.obs["coarse"]
qry.obs["coarse"] = "Unknown"  # surgery sees no query labels

# scArches-friendly architecture (Lotfollahi et al. 2022)
arch = dict(n_layers=2, n_latent=30, encode_covariates=True, deeply_inject_covariates=False,
            use_layer_norm="both", use_batch_norm="none")
scvi.model.SCVI.setup_anndata(ref, layer="counts", batch_key=BATCH)
vae = scvi.model.SCVI(ref, **arch)
vae.train(max_epochs=200, early_stopping=True)
lvae = scvi.model.SCANVI.from_scvi_model(vae, labels_key="coarse", unlabeled_category="Unknown")
lvae.train(max_epochs=20, n_samples_per_label=100)
for m in (vae, lvae):
    m.save(f"models/{type(m).__name__.lower()}_ref", overwrite=True)

surgery = dict(max_epochs=100, plan_kwargs=dict(weight_decay=0.0))
q_vae = scvi.model.SCVI.load_query_data(qry, vae)
q_vae.train(**surgery)
q_lvae = scvi.model.SCANVI.load_query_data(qry, lvae)
q_lvae.train(**surgery)

for name, m, qm in (("scvi", vae, q_vae), ("scanvi", lvae, q_lvae)):
    emb = np.vstack([m.get_latent_representation(ref), qm.get_latent_representation(qry)])
    np.save(f"results/emb_{name}.npy", emb.astype(np.float32))

soft = q_lvae.predict(qry, soft=True)
pd.DataFrame({"pred": soft.idxmax(axis=1), "conf": soft.max(axis=1)}, index=qry.obs_names).to_csv("results/pred_scanvi.csv")
