import pymde
import numpy as np
import pandas as pd
import plotnine as gg
import patchworklib as pw
pymde.seed(0)

import scanpy as sc
from scvi.model import SCVI, SCANVI

from typing import Tuple

def immediate_early_pc(adata) -> np.ndarray:
    """
    Returns the first principal component of various immediate early genes.
    """
    adata = adata.copy()  
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    ieg_list = ["Fos","Fosb","Jun","Junb","Egr1","Egr2","Npas4","Nr4a1","Nr4a2","Nr4a3"]
    adata = adata[:, ieg_list].copy()
    sc.pp.scale(adata)

    ieg_pc = sc.pp.pca(adata.X, n_comps=1)
    ieg_pc = np.array(ieg_pc).flatten()

    return ieg_pc

def integrate_datasets(adata, n_genes, early_stopping = True) -> Tuple[SCVI, SCANVI]:
    adata = adata.copy()

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes = n_genes, subset=True, batch_key="Dataset")

    setup_kwargs = {
        "layer": "counts",
        "batch_key": "Sample_ID",
        "categorical_covariate_keys": [
            "suspension_type",
            "assay",
            "sex",
        ],
        "continuous_covariate_keys": [
            "immediate_early",
            "pctMito_RNA",
        ],
    }
    scvi_init = {
        "n_layers": 3,
        "n_latent": 20,
        "dropout_rate": 0.01,
        "dispersion": "gene-cell",
        "gene_likelihood": "zinb",
    }
    scanvi_init = {
        "labels_key": "C66_named",
        "unlabeled_category": "nan",
    }
    plan_kwargs = {
        "lr": 0.01,
        "lr_min": 0.0001,
        "reduce_lr_on_plateau": True,
        "lr_patience": 8,
        "lr_factor": 0.1**0.25,
        "lr_scheduler_metric": "elbo_validation",
    }
    train_kwargs = {
        "check_val_every_n_epoch": 1,
        "early_stopping": early_stopping,
        "early_stopping_patience": 25,
        "early_stopping_monitor": "elbo_validation",
    }
    scvi_train = {
        "max_epochs": 400,
        "plan_kwargs": {
            "n_epochs_kl_warmup": 400,
            **plan_kwargs
        },
        **train_kwargs,
    }
    scanvi_train = {
        "max_epochs": 250,
        "plan_kwargs": {
            "classification_ratio": 200,
            "n_epochs_kl_warmup": 100,
            **plan_kwargs,
        },
        **train_kwargs,
    }

    label = adata.obs[scanvi_init["labels_key"]]
    label = np.array(label, dtype = str)
    adata.obs[scanvi_init["labels_key"]] = label

    SCVI.setup_anndata(adata, **setup_kwargs) 

    scvae = SCVI(adata, **scvi_init)
    scvae.train(**scvi_train)
    
    scanvae = SCANVI.from_scvi_model(scvae, **scanvi_init)
    scanvae.train(**scanvi_train)

    return scvae, scanvae

def plot_loss(vae):
    def _validation_len() -> int:
        return vae.validation_indices.__len__()

    def _as_array(key: str) -> np.ndarray:
        return np.array(vae.history[key], dtype = float).flatten()
    
    df = pd.DataFrame({
        "index": np.arange(len(_as_array("validation_loss"))),
        "validation": _as_array("validation_loss"),
        "elbo": _as_array("elbo_validation"),
        "reconstruction": _as_array("reconstruction_loss_validation") / _validation_len(),
        "kldiv": _as_array("kl_local_validation") / _validation_len(),
        **{
            "classification": _as_array("classification_loss_validation") 
            if "classification_loss_validation" in vae.history.keys() 
            else np.zeros_like(_as_array("elbo_validation")),
        }
    })
    df["custom"] = df["elbo"] + df["classification"] * 200

    def _plot_loss(df, key: str):
        p = gg.ggplot(df, gg.aes(x = "index")) + \
            gg.geom_line(gg.aes(y = key), color = "black") + \
            gg.scale_x_log10() + \
            gg.theme_bw()
        return pw.load_ggplot(p, figsize=(3,3))

    p = (_plot_loss(df,"validation")|_plot_loss(df, "elbo"))/(_plot_loss(df, "classification")|_plot_loss(df, "custom"))

    return p

def process_integration(adata: sc.AnnData, vae: SCANVI) -> None:
    # Get SCANVI predictions
    adata.obs["C66_predicted"] = vae.predict()
    
    # Get latent representation
    vae = vae.get_latent_representation()
    adata.obsm["X_vae"] = vae


def mde_embedding(adata: sc.AnnData) -> None:  
    vae = adata.obsm["X_vae"]
    mde_kwargs = {
        'attractive_penalty': pymde.penalties.Huber,
        'repulsive_penalty': pymde.penalties.Log,
        'n_neighbors': 30, 
        'device': 'cuda'
    }
    mde = pymde.preserve_neighbors(vae, **mde_kwargs)
    mde = np.array(mde.embed().cpu())

    adata.obsm["X_mde"] = mde

