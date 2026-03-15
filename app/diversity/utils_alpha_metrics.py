"""Alpha diversity metric calculations."""

import numpy as np
from scipy import stats as sp_stats


def observed_features(counts):
    """Count of features with count > 0."""
    return int(np.sum(np.asarray(counts) > 0))


def chao1(counts):
    """Chao1 estimated richness."""
    counts = np.asarray(counts, dtype=float)
    s_obs = np.sum(counts > 0)
    n1 = np.sum(counts == 1)  # singletons
    n2 = np.sum(counts == 2)  # doubletons

    if n2 == 0:
        # Bias-corrected form
        return float(s_obs + n1 * (n1 - 1) / 2) if n1 > 0 else float(s_obs)

    return float(s_obs + (n1 ** 2) / (2 * n2))


def ace(counts):
    """Abundance-based Coverage Estimator."""
    counts = np.asarray(counts, dtype=float)
    counts = counts[counts > 0]

    # Rare species (count <= 10) and abundant species
    rare = counts[counts <= 10]
    s_rare = len(rare)
    s_abund = len(counts[counts > 10])

    if s_rare == 0:
        return float(s_abund)

    n_rare = rare.sum()
    if n_rare == 0:
        return float(s_abund)

    n1 = np.sum(rare == 1)
    c_ace = 1 - n1 / n_rare

    if c_ace == 0:
        return float(s_abund + s_rare)

    # Coefficient of variation
    sum_fi_i = np.sum(rare * (rare - 1))
    gamma2 = max((s_rare * sum_fi_i) / (c_ace * n_rare * (n_rare - 1)) - 1, 0)

    return float(s_abund + s_rare / c_ace + n1 * gamma2 / c_ace)


def shannon(counts):
    """Shannon diversity index H'."""
    counts = np.asarray(counts, dtype=float)
    counts = counts[counts > 0]
    total = counts.sum()
    if total == 0:
        return 0.0
    p = counts / total
    return float(-np.sum(p * np.log(p)))


def pielou_evenness(counts):
    """Pielou's evenness J = H' / ln(S)."""
    h = shannon(counts)
    s = observed_features(counts)
    if s <= 1:
        return 0.0
    return float(h / np.log(s))


def simpson(counts):
    """Simpson index (1 - D)."""
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total <= 1:
        return 0.0
    p = counts / total
    return float(1 - np.sum(p ** 2))


def simpson_evenness(counts):
    """Simpson evenness = 1 / (D * S)."""
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total <= 1:
        return 0.0
    p = counts / total
    d = np.sum(p ** 2)
    s = np.sum(counts > 0)
    if d == 0 or s == 0:
        return 0.0
    return float(1 / (d * s))


def inverse_simpson(counts):
    """Inverse Simpson 1/D."""
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total <= 1:
        return 0.0
    p = counts / total
    d = np.sum(p ** 2)
    if d == 0:
        return 0.0
    return float(1 / d)


def berger_parker(counts):
    """Berger-Parker dominance: max(n_i) / N."""
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()
    if total == 0:
        return 0.0
    return float(counts.max() / total)


def faiths_pd(counts, tree, feature_ids):
    """
    Faith's Phylogenetic Diversity.
    Requires a scikit-bio TreeNode and matching feature IDs.
    """
    from skbio.diversity import alpha_diversity

    counts = np.asarray(counts, dtype=int)
    result = alpha_diversity(
        "faith_pd", [counts], ids=[feature_ids], tree=tree, otu_ids=feature_ids
    )
    return float(result.iloc[0])


def calc_all_metrics(abundance_df, tree=None):
    """
    Calculate all alpha diversity metrics for all samples.
    Returns dict of metric_name -> list of values (one per sample).
    """
    samples = list(abundance_df.columns)
    feature_ids = list(abundance_df.index.astype(str))

    metrics = {
        "observed_features": [],
        "chao1": [],
        "ace": [],
        "shannon": [],
        "pielou_evenness": [],
        "simpson": [],
        "simpson_evenness": [],
        "inverse_simpson": [],
        "berger_parker": [],
    }

    for sample in samples:
        counts = abundance_df[sample].values

        metrics["observed_features"].append(observed_features(counts))
        metrics["chao1"].append(round(chao1(counts), 2))
        metrics["ace"].append(round(ace(counts), 2))
        metrics["shannon"].append(round(shannon(counts), 4))
        metrics["pielou_evenness"].append(round(pielou_evenness(counts), 4))
        metrics["simpson"].append(round(simpson(counts), 4))
        metrics["simpson_evenness"].append(round(simpson_evenness(counts), 4))
        metrics["inverse_simpson"].append(round(inverse_simpson(counts), 2))
        metrics["berger_parker"].append(round(berger_parker(counts), 4))

    # Faith's PD if tree provided
    if tree is not None:
        try:
            from skbio.diversity import alpha_diversity as skbio_alpha

            counts_matrix = abundance_df.T.values.astype(int)
            result = skbio_alpha(
                "faith_pd", counts_matrix,
                ids=samples, tree=tree, otu_ids=feature_ids,
            )
            metrics["faiths_pd"] = [round(float(v), 2) for v in result.values]
        except Exception:
            metrics["faiths_pd"] = None
    else:
        metrics["faiths_pd"] = None

    return metrics


METRIC_LABELS = {
    "observed_features": "Observed Features",
    "chao1": "Chao1",
    "ace": "ACE",
    "shannon": "Shannon (H')",
    "pielou_evenness": "Pielou's Evenness (J)",
    "simpson": "Simpson (1-D)",
    "simpson_evenness": "Simpson Evenness",
    "inverse_simpson": "Inverse Simpson",
    "berger_parker": "Berger-Parker",
    "faiths_pd": "Faith's PD",
}
