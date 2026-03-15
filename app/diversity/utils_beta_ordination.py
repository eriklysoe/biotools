"""Beta diversity ordination: PCoA, NMDS, t-SNE, PCA."""

import numpy as np


def run_pcoa(dm):
    """Run PCoA ordination. Returns dict with coordinates and variance explained."""
    from skbio.stats.ordination import pcoa

    result = pcoa(dm)
    samples = list(dm.ids)

    coords = result.samples
    proportion = result.proportion_explained

    n_axes = min(3, coords.shape[1])

    return {
        "method": "pcoa",
        "axis1": coords.iloc[:, 0].tolist(),
        "axis2": coords.iloc[:, 1].tolist() if n_axes >= 2 else [0] * len(samples),
        "axis3": coords.iloc[:, 2].tolist() if n_axes >= 3 else [0] * len(samples),
        "variance_explained": [
            round(float(proportion.iloc[i]) * 100, 1)
            for i in range(n_axes)
        ],
        "stress": None,
    }


def run_nmds(dm, n_init=100):
    """Run NMDS. Returns dict with coordinates and stress."""
    from sklearn.manifold import MDS

    matrix = dm.data

    nmds = MDS(
        n_components=2,
        metric=False,
        dissimilarity="precomputed",
        n_init=n_init,
        max_iter=300,
        random_state=42,
    )

    coords = nmds.fit_transform(matrix)

    return {
        "method": "nmds",
        "axis1": coords[:, 0].tolist(),
        "axis2": coords[:, 1].tolist(),
        "axis3": [0.0] * len(dm.ids),
        "variance_explained": None,
        "stress": round(float(nmds.stress_), 4),
    }


def run_tsne(dm, perplexity=30):
    """Run t-SNE. Returns dict with coordinates."""
    from sklearn.manifold import TSNE

    matrix = dm.data
    perplexity = min(perplexity, max(1, len(dm.ids) - 1))

    tsne = TSNE(
        n_components=2,
        metric="precomputed",
        perplexity=perplexity,
        random_state=42,
        init="random",
    )

    coords = tsne.fit_transform(matrix)

    return {
        "method": "tsne",
        "axis1": coords[:, 0].tolist(),
        "axis2": coords[:, 1].tolist(),
        "axis3": [0.0] * len(dm.ids),
        "variance_explained": None,
        "stress": None,
    }


def run_pca_clr(abundance_df):
    """
    Run PCA on CLR-transformed data (Aitchison PCA).
    Returns dict with coordinates, variance explained, and loadings.
    """
    from sklearn.decomposition import PCA

    # CLR transform
    data = abundance_df.values.astype(float) + 1
    log_data = np.log(data)
    geometric_mean = log_data.mean(axis=0)
    clr = log_data - geometric_mean

    # PCA on transposed data (samples x features)
    clr_t = clr.T

    n_components = min(3, clr_t.shape[0], clr_t.shape[1])
    pca = PCA(n_components=n_components, random_state=42)
    coords = pca.fit_transform(clr_t)

    # Feature loadings for biplot
    loadings = None
    feature_ids = list(abundance_df.index.astype(str))
    if len(feature_ids) <= 50:
        loadings = {
            "feature_ids": feature_ids,
            "pc1": pca.components_[0].tolist(),
            "pc2": pca.components_[1].tolist() if n_components >= 2 else [0] * len(feature_ids),
        }

    return {
        "method": "pca_clr",
        "axis1": coords[:, 0].tolist(),
        "axis2": coords[:, 1].tolist() if n_components >= 2 else [0] * len(abundance_df.columns),
        "axis3": coords[:, 2].tolist() if n_components >= 3 else [0] * len(abundance_df.columns),
        "variance_explained": [
            round(float(v) * 100, 1) for v in pca.explained_variance_ratio_[:n_components]
        ],
        "stress": None,
        "loadings": loadings,
    }


def run_ordination(abundance_df, dm, method, perplexity=30):
    """Dispatch to the correct ordination method."""
    if method == "pcoa":
        return run_pcoa(dm)
    elif method == "nmds":
        return run_nmds(dm)
    elif method == "tsne":
        return run_tsne(dm, perplexity=perplexity)
    elif method == "pca_clr":
        return run_pca_clr(abundance_df)
    else:
        raise ValueError(f"Unknown ordination method: {method}")
