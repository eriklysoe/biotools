"""Beta diversity distance matrix calculation."""

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform


DISTANCE_METRICS = {
    "compositional": {
        "label": "Compositional",
        "metrics": {
            "braycurtis": "Bray-Curtis",
            "jaccard": "Jaccard (presence/absence)",
            "canberra": "Canberra",
        },
    },
    "phylogenetic": {
        "label": "Phylogenetic (requires tree)",
        "metrics": {
            "unweighted_unifrac": "Unweighted UniFrac",
            "weighted_unifrac": "Weighted UniFrac",
        },
    },
    "general": {
        "label": "General",
        "metrics": {
            "euclidean": "Euclidean",
            "correlation": "Correlation (1 - Pearson)",
            "aitchison": "Aitchison (CLR + Euclidean)",
        },
    },
}


def get_flat_metrics():
    """Return flat dict of metric_id -> label."""
    flat = {}
    for group in DISTANCE_METRICS.values():
        flat.update(group["metrics"])
    return flat


def _clr_transform(df):
    """Centered log-ratio transformation. Adds pseudocount of 1."""
    data = df.values.astype(float) + 1  # pseudocount
    log_data = np.log(data)
    geometric_mean = log_data.mean(axis=0)
    clr = log_data - geometric_mean
    return pd.DataFrame(clr, index=df.index, columns=df.columns)


def compute_distance_matrix(abundance_df, metric, tree=None):
    """
    Compute distance matrix for given metric.

    Args:
        abundance_df: features x samples DataFrame
        metric: distance metric string
        tree: optional skbio TreeNode

    Returns:
        (skbio DistanceMatrix, warnings list)
    """
    from skbio import DistanceMatrix

    warnings = []
    samples = list(abundance_df.columns.astype(str))

    # Transpose: samples x features for distance calculations
    data = abundance_df.T.values.astype(float)

    if metric in ("unweighted_unifrac", "weighted_unifrac"):
        if tree is None:
            raise ValueError(f"{metric} requires a phylogenetic tree")

        from skbio.diversity import beta_diversity

        feature_ids = list(abundance_df.index.astype(str))
        counts = abundance_df.T.values.astype(int)

        skbio_metric = metric.replace("_", "-")
        dm = beta_diversity(
            skbio_metric, counts,
            ids=samples, otu_ids=feature_ids, tree=tree,
        )
        return dm, warnings

    elif metric == "aitchison":
        # CLR transform then Euclidean
        n_zeros = int((abundance_df.values == 0).sum())
        if n_zeros > 0:
            warnings.append(
                f"Aitchison distance: pseudocount of 1 added to {n_zeros} zero values"
            )
        clr_df = _clr_transform(abundance_df)
        data = clr_df.T.values
        dists = pdist(data, metric="euclidean")

    elif metric == "jaccard":
        # Binary Jaccard
        binary = (data > 0).astype(float)
        dists = pdist(binary, metric="jaccard")

    elif metric in ("braycurtis", "canberra", "euclidean", "correlation"):
        dists = pdist(data, metric=metric)

    else:
        raise ValueError(f"Unknown distance metric: {metric}")

    # Handle NaN
    if np.any(np.isnan(dists)):
        warnings.append("Distance matrix contains NaN values — some samples may be identical or have zero variance")
        dists = np.nan_to_num(dists, nan=0.0)

    matrix = squareform(dists)
    dm = DistanceMatrix(matrix, ids=samples)

    return dm, warnings
