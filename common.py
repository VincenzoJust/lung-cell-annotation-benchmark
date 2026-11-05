import anndata as ad
import scanpy as sc

BATCH = "dataset"  # HLCA study; query cells get their donor_id here


def load():
    """Raw-count ref and query, with log1p-normalized X and counts in layers['counts']."""
    ref, qry = ad.read_h5ad("data/ref.h5ad"), ad.read_h5ad("data/query.h5ad")
    qry.obs[BATCH] = qry.obs["donor_id"].astype(str)
    for a in (ref, qry):
        a.obs[BATCH] = a.obs[BATCH].astype(str)
        a.layers["counts"] = a.X.copy()
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    return ref, qry


def hvg(ref, n=2000):
    sc.pp.highly_variable_genes(ref, n_top_genes=n, flavor="seurat_v3", layer="counts")
    return ref.var_names[ref.var.highly_variable]
