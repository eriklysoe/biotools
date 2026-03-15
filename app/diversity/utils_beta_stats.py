"""Beta diversity statistical testing: PERMANOVA, PERMDISP, ANOSIM, Mantel."""

import numpy as np
import pandas as pd


def _benjamini_hochberg(p_values):
    """Benjamini-Hochberg FDR correction."""
    n = len(p_values)
    if n == 0:
        return []
    p_arr = np.array(p_values)
    ranked = np.argsort(p_arr)
    adjusted = np.zeros(n)
    for i, idx in enumerate(reversed(ranked)):
        rank = n - i
        if i == 0:
            adjusted[idx] = min(p_arr[idx] * n / rank, 1.0)
        else:
            adjusted[idx] = min(p_arr[idx] * n / rank, adjusted[ranked[n - i]])
    return adjusted.tolist()


def run_permanova(dm, grouping, permutations=999):
    """
    Run PERMANOVA test.
    Returns dict with F-statistic, R², p-value.
    """
    from skbio.stats.distance import permanova

    result = permanova(dm, grouping, permutations=permutations)

    return {
        "test": "PERMANOVA",
        "f_statistic": round(float(result["test statistic"]), 4),
        "p_value": round(float(result["p-value"]), 6),
        "permutations": int(result["number of permutations"]),
        "sample_size": int(result["sample size"]),
    }


def run_permdisp(dm, grouping, permutations=999):
    """
    Run PERMDISP (betadisper) test.
    Returns dict with F-statistic, p-value, warning if significant.
    """
    from skbio.stats.distance import permdisp

    result = permdisp(dm, grouping, permutations=permutations)

    p = float(result["p-value"])
    warning = None
    if p < 0.05:
        warning = (
            "PERMDISP is significant (p={:.4f}), indicating unequal group dispersions. "
            "PERMANOVA result may reflect dispersion differences rather than "
            "centroid differences.".format(p)
        )

    return {
        "test": "PERMDISP",
        "f_statistic": round(float(result["test statistic"]), 4),
        "p_value": round(p, 6),
        "permutations": int(result["number of permutations"]),
        "warning": warning,
    }


def run_anosim(dm, grouping, permutations=999):
    """
    Run ANOSIM test.
    Returns dict with R statistic, p-value.
    """
    from skbio.stats.distance import anosim

    result = anosim(dm, grouping, permutations=permutations)

    return {
        "test": "ANOSIM",
        "r_statistic": round(float(result["test statistic"]), 4),
        "p_value": round(float(result["p-value"]), 6),
        "permutations": int(result["number of permutations"]),
    }


def run_pairwise_permanova(dm, grouping, permutations=999):
    """
    Run pairwise PERMANOVA between all group pairs.
    Returns list of comparison dicts with BH-corrected p-values.
    """
    from skbio.stats.distance import permanova

    groups = sorted(set(grouping))
    if len(groups) < 3:
        return None

    samples = list(dm.ids)
    comparisons = []
    p_values = []

    for i, g1 in enumerate(groups):
        for g2 in groups[i+1:]:
            # Subset to just these two groups
            mask = [(grouping[j] == g1 or grouping[j] == g2) for j in range(len(grouping))]
            subset_ids = [samples[j] for j in range(len(samples)) if mask[j]]
            subset_grouping = [grouping[j] for j in range(len(grouping)) if mask[j]]

            if len(subset_ids) < 3:
                continue

            sub_dm = dm.filter(subset_ids)
            try:
                result = permanova(sub_dm, subset_grouping, permutations=permutations)
                f_stat = float(result["test statistic"])
                p_val = float(result["p-value"])
            except Exception:
                f_stat = 0.0
                p_val = 1.0

            p_values.append(p_val)
            comparisons.append({
                "group1": str(g1),
                "group2": str(g2),
                "f_statistic": round(f_stat, 4),
                "p_value": round(p_val, 6),
            })

    # BH correction
    adjusted = _benjamini_hochberg(p_values)
    for comp, p_adj in zip(comparisons, adjusted):
        comp["p_adjusted"] = round(p_adj, 6)
        comp["significant"] = p_adj < 0.05

    return comparisons


def run_mantel(dm1, dm2, permutations=999):
    """
    Run Mantel test between two distance matrices.
    Returns dict with r statistic, p-value.
    """
    from skbio.stats.distance import mantel

    r, p, n = mantel(dm1, dm2, permutations=permutations)

    return {
        "test": "Mantel",
        "r_statistic": round(float(r), 4),
        "p_value": round(float(p), 6),
        "permutations": permutations,
        "n": int(n),
    }


def run_all_beta_stats(dm, grouping_series, permutations=999):
    """
    Run all beta diversity stats: PERMANOVA, PERMDISP, ANOSIM, pairwise PERMANOVA.
    grouping_series: list of group labels aligned with dm.ids.
    """
    groups = sorted(set(grouping_series))
    warnings = []

    # Check group sizes
    from collections import Counter
    sizes = Counter(grouping_series)
    min_size = min(sizes.values())
    max_size = max(sizes.values())
    if max_size > 3 * min_size:
        warnings.append(
            f"Unequal group sizes (range: {min_size}-{max_size}) may affect PERMANOVA results"
        )

    result = {
        "grouping_column": None,  # set by caller
        "warnings": warnings,
    }

    try:
        result["permanova"] = run_permanova(dm, grouping_series, permutations)
    except Exception as e:
        result["permanova"] = {"error": str(e)}

    try:
        result["permdisp"] = run_permdisp(dm, grouping_series, permutations)
    except Exception as e:
        result["permdisp"] = {"error": str(e)}

    try:
        result["anosim"] = run_anosim(dm, grouping_series, permutations)
    except Exception as e:
        result["anosim"] = {"error": str(e)}

    if len(groups) >= 3:
        try:
            result["pairwise_permanova"] = run_pairwise_permanova(
                dm, grouping_series, permutations
            )
        except Exception as e:
            result["pairwise_permanova"] = {"error": str(e)}
    else:
        result["pairwise_permanova"] = None

    return result
