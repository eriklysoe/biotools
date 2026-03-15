"""Alpha diversity statistical testing."""

import numpy as np
from scipy import stats as sp_stats


def _significance_stars(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    return "ns"


def _benjamini_hochberg(p_values):
    """Benjamini-Hochberg FDR correction."""
    n = len(p_values)
    if n == 0:
        return []
    ranked = np.argsort(p_values)
    adjusted = np.zeros(n)
    for i, idx in enumerate(reversed(ranked)):
        rank = n - i
        if i == 0:
            adjusted[idx] = min(p_values[idx] * n / rank, 1.0)
        else:
            adjusted[idx] = min(p_values[idx] * n / rank, adjusted[ranked[n - i]])
    return adjusted.tolist()


def run_alpha_stats(values, groups, paired=False):
    """
    Run statistical tests on alpha diversity values by group.

    Args:
        values: list of diversity values (one per sample)
        groups: list of group labels (one per sample)

    Returns:
        dict with test results
    """
    values = np.array(values, dtype=float)
    groups = np.array(groups)

    unique_groups = sorted(set(groups))
    n_groups = len(unique_groups)

    if n_groups < 2:
        return {"error": "Need at least 2 groups for statistical testing"}

    group_data = {g: values[groups == g] for g in unique_groups}

    result = {}

    if n_groups == 2:
        g1, g2 = unique_groups
        d1, d2 = group_data[g1], group_data[g2]

        if paired and len(d1) == len(d2):
            stat, p = sp_stats.wilcoxon(d1, d2)
            result["test"] = "wilcoxon"
        else:
            stat, p = sp_stats.mannwhitneyu(d1, d2, alternative="two-sided")
            result["test"] = "mann-whitney"

        result["statistic"] = round(float(stat), 4)
        result["p_value"] = round(float(p), 6)
        result["significance"] = _significance_stars(p)
        result["posthoc"] = None

    else:
        # Kruskal-Wallis
        group_arrays = [group_data[g] for g in unique_groups]
        stat, p = sp_stats.kruskal(*group_arrays)

        result["test"] = "kruskal-wallis"
        result["statistic"] = round(float(stat), 4)
        result["p_value"] = round(float(p), 6)
        result["significance"] = _significance_stars(p)

        # Dunn's post-hoc if significant
        if p < 0.05:
            result["posthoc"] = _dunns_test(values, groups, unique_groups)
        else:
            result["posthoc"] = None

    # Group summaries
    result["group_summary"] = {}
    for g in unique_groups:
        d = group_data[g]
        result["group_summary"][g] = {
            "n": int(len(d)),
            "mean": round(float(np.mean(d)), 4),
            "median": round(float(np.median(d)), 4),
            "sd": round(float(np.std(d, ddof=1)), 4) if len(d) > 1 else 0,
            "min": round(float(np.min(d)), 4),
            "max": round(float(np.max(d)), 4),
        }

    return result


def _dunns_test(values, groups, unique_groups):
    """Dunn's post-hoc test with BH correction."""
    try:
        import scikit_posthocs as sp
        import pandas as pd

        result = sp.posthoc_dunn(
            pd.DataFrame({"value": values, "group": groups}),
            val_col="value", group_col="group", p_adjust="fdr_bh",
        )

        comparisons = []
        for i, g1 in enumerate(unique_groups):
            for g2 in unique_groups[i+1:]:
                p_adj = float(result.loc[g1, g2])
                comparisons.append({
                    "group1": str(g1),
                    "group2": str(g2),
                    "p_adjusted": round(p_adj, 6),
                    "significant": p_adj < 0.05,
                    "significance": _significance_stars(p_adj),
                })

        return comparisons

    except ImportError:
        # Fallback: manual pairwise Mann-Whitney with BH correction
        comparisons = []
        p_values = []

        for i, g1 in enumerate(unique_groups):
            for g2 in unique_groups[i+1:]:
                d1 = values[groups == g1]
                d2 = values[groups == g2]
                _, p = sp_stats.mannwhitneyu(d1, d2, alternative="two-sided")
                p_values.append(p)
                comparisons.append({
                    "group1": str(g1),
                    "group2": str(g2),
                    "p_value": round(float(p), 6),
                })

        adjusted = _benjamini_hochberg(np.array(p_values))
        for comp, p_adj in zip(comparisons, adjusted):
            comp["p_adjusted"] = round(p_adj, 6)
            comp["significant"] = p_adj < 0.05
            comp["significance"] = _significance_stars(p_adj)

        return comparisons
