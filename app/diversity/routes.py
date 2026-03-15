"""Diversity routes – Alpha & Beta diversity endpoints."""

import os
import traceback

from flask import jsonify, render_template, request

from . import diversity_bp
from .utils_shared import (
    parse_abundance_table, parse_metadata, match_metadata,
    get_sample_depths, rarefy_table, parse_newick_tree, validate_tree_tips,
)
from .utils_alpha_metrics import calc_all_metrics, METRIC_LABELS
from .utils_alpha_stats import run_alpha_stats
from .utils_alpha_rarefaction import calc_rarefaction_curves
from .utils_beta_distances import compute_distance_matrix, DISTANCE_METRICS, get_flat_metrics
from .utils_beta_ordination import run_ordination
from .utils_beta_stats import run_all_beta_stats, run_mantel

BASE_URL = os.environ.get("BASE_URL", "http://localhost:5590").rstrip("/")


@diversity_bp.route("/diversity")
def diversity_page():
    return render_template(
        "diversity.html",
        base_url=BASE_URL,
        distance_metrics=DISTANCE_METRICS,
        metric_labels=METRIC_LABELS,
    )


# ═══════ Alpha Diversity ═══════

@diversity_bp.route("/api/diversity/alpha", methods=["POST"])
def api_alpha():
    """Calculate alpha diversity metrics."""
    if "abundance_file" not in request.files or not request.files["abundance_file"].filename:
        return jsonify({"error": "No abundance table uploaded"}), 400

    abundance_bytes = request.files["abundance_file"].read()
    all_warnings = []

    # Parse abundance table
    try:
        abundance_df, parse_warnings = parse_abundance_table(abundance_bytes)
        all_warnings.extend(parse_warnings)
    except ValueError as e:
        return jsonify({"error": str(e)}), 400

    # Parse optional metadata
    metadata_df = None
    metadata_columns = []
    metadata_dict = {}
    if "metadata_file" in request.files and request.files["metadata_file"].filename:
        try:
            meta_bytes = request.files["metadata_file"].read()
            metadata_df = parse_metadata(meta_bytes)
            metadata_df, meta_warnings = match_metadata(abundance_df, metadata_df)
            all_warnings.extend(meta_warnings)
            metadata_columns = list(metadata_df.columns)
            # Align metadata to abundance columns
            for col in metadata_columns:
                metadata_dict[col] = [
                    str(metadata_df.loc[s, col]) if s in metadata_df.index else "NA"
                    for s in abundance_df.columns.astype(str)
                ]
        except ValueError as e:
            all_warnings.append(f"Metadata error: {e}")

    # Parse optional tree
    tree = None
    if "tree_file" in request.files and request.files["tree_file"].filename:
        try:
            tree_bytes = request.files["tree_file"].read()
            tree = parse_newick_tree(tree_bytes)
            ok, msg = validate_tree_tips(tree, abundance_df.index)
            if not ok:
                all_warnings.append(msg)
                tree = None
            elif msg:
                all_warnings.append(msg)
        except Exception as e:
            all_warnings.append(f"Tree error: {e}")
            tree = None

    # Rarefaction
    seed = 42
    rarefaction_depth = request.form.get("rarefaction_depth")
    sample_depths = get_sample_depths(abundance_df)
    rarefied = False

    if rarefaction_depth:
        try:
            rarefaction_depth = int(rarefaction_depth)
            abundance_df, dropped, rare_warnings = rarefy_table(
                abundance_df, rarefaction_depth, seed=seed
            )
            all_warnings.extend(rare_warnings)
            rarefied = True
            # Update metadata alignment if samples were dropped
            if dropped and metadata_dict:
                remaining = list(abundance_df.columns.astype(str))
                for col in metadata_columns:
                    metadata_dict[col] = [
                        str(metadata_df.loc[s, col]) if s in metadata_df.index else "NA"
                        for s in remaining
                    ]
        except ValueError as e:
            return jsonify({"error": str(e)}), 400
    else:
        rarefaction_depth = None

    # Calculate metrics
    try:
        metrics = calc_all_metrics(abundance_df, tree=tree)
    except Exception:
        traceback.print_exc()
        return jsonify({"error": "Error calculating diversity metrics"}), 500

    # Statistical testing
    statistics = {}
    grouping_col = request.form.get("grouping_column")
    paired = request.form.get("paired", "false").lower() == "true"

    if grouping_col and grouping_col in metadata_dict:
        groups = metadata_dict[grouping_col]
        unique_groups = set(groups)
        if len(unique_groups) >= 2:
            for metric_name, values in metrics.items():
                if values is None:
                    continue
                try:
                    stat_result = run_alpha_stats(values, groups, paired=paired)
                    statistics[metric_name] = stat_result
                except Exception:
                    pass

    # Rarefaction curves (on original un-rarefied data)
    try:
        abundance_bytes_for_rare = abundance_bytes
        rare_df, _ = parse_abundance_table(abundance_bytes_for_rare)
        rarefaction_curves = calc_rarefaction_curves(rare_df, seed=seed)
    except Exception:
        rarefaction_curves = None

    samples = [str(s) for s in abundance_df.columns]

    return jsonify({
        "samples": samples,
        "metrics": metrics,
        "metric_labels": METRIC_LABELS,
        "metadata": metadata_dict if metadata_dict else None,
        "metadata_columns": metadata_columns,
        "rarefaction_depth": rarefaction_depth,
        "rarefaction_seed": seed,
        "sample_depths": sample_depths,
        "statistics": statistics if statistics else None,
        "rarefaction_curves": rarefaction_curves,
        "num_samples": len(samples),
        "num_features": int(abundance_df.shape[0]),
        "warnings": all_warnings,
    })


# ═══════ Beta Diversity ═══════

@diversity_bp.route("/api/diversity/beta", methods=["POST"])
def api_beta():
    """Calculate beta diversity with ordination and statistics."""
    if "abundance_file" not in request.files or not request.files["abundance_file"].filename:
        return jsonify({"error": "No abundance table uploaded"}), 400

    abundance_bytes = request.files["abundance_file"].read()
    all_warnings = []

    # Parse abundance table
    try:
        abundance_df, parse_warnings = parse_abundance_table(abundance_bytes)
        all_warnings.extend(parse_warnings)
    except ValueError as e:
        return jsonify({"error": str(e)}), 400

    if abundance_df.shape[1] < 3:
        return jsonify({"error": "Need at least 3 samples for beta diversity"}), 400

    # Parse optional metadata
    metadata_df = None
    metadata_columns = []
    metadata_dict = {}
    if "metadata_file" in request.files and request.files["metadata_file"].filename:
        try:
            meta_bytes = request.files["metadata_file"].read()
            metadata_df = parse_metadata(meta_bytes)
            metadata_df, meta_warnings = match_metadata(abundance_df, metadata_df)
            all_warnings.extend(meta_warnings)
            metadata_columns = list(metadata_df.columns)
        except ValueError as e:
            all_warnings.append(f"Metadata error: {e}")

    # Parse optional tree
    tree = None
    if "tree_file" in request.files and request.files["tree_file"].filename:
        try:
            tree_bytes = request.files["tree_file"].read()
            tree = parse_newick_tree(tree_bytes)
            ok, msg = validate_tree_tips(tree, abundance_df.index)
            if not ok:
                all_warnings.append(msg)
                tree = None
            elif msg:
                all_warnings.append(msg)
        except Exception as e:
            all_warnings.append(f"Tree error: {e}")
            tree = None

    # Rarefaction
    seed = 42
    rarefaction_depth = request.form.get("rarefaction_depth")
    sample_depths = get_sample_depths(abundance_df)

    if rarefaction_depth:
        try:
            rarefaction_depth = int(rarefaction_depth)
            abundance_df, dropped, rare_warnings = rarefy_table(
                abundance_df, rarefaction_depth, seed=seed
            )
            all_warnings.extend(rare_warnings)
        except ValueError as e:
            return jsonify({"error": str(e)}), 400
    else:
        rarefaction_depth = None

    samples = [str(s) for s in abundance_df.columns]

    # Align metadata
    if metadata_df is not None:
        for col in metadata_columns:
            metadata_dict[col] = [
                str(metadata_df.loc[s, col]) if s in metadata_df.index else "NA"
                for s in samples
            ]

    # Distance metric
    metric = request.form.get("distance_metric", "braycurtis")
    if metric in ("unweighted_unifrac", "weighted_unifrac") and tree is None:
        return jsonify({"error": f"{metric} requires a phylogenetic tree file"}), 400

    try:
        dm, dist_warnings = compute_distance_matrix(abundance_df, metric, tree=tree)
        all_warnings.extend(dist_warnings)
    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    except Exception:
        traceback.print_exc()
        return jsonify({"error": "Error computing distance matrix"}), 500

    # Ordination
    ord_method = request.form.get("ordination_method", "pcoa")
    perplexity = int(request.form.get("perplexity", 30))

    try:
        ordination = run_ordination(abundance_df, dm, ord_method, perplexity=perplexity)
    except Exception:
        traceback.print_exc()
        return jsonify({"error": f"Ordination failed ({ord_method})"}), 500

    # Statistics
    statistics = None
    grouping_col = request.form.get("grouping_column")
    permutations = int(request.form.get("permutations", 999))

    if grouping_col and grouping_col in metadata_dict:
        grouping_series = metadata_dict[grouping_col]
        unique_groups = set(grouping_series)
        if len(unique_groups) >= 2:
            try:
                statistics = run_all_beta_stats(dm, grouping_series, permutations)
                statistics["grouping_column"] = grouping_col
            except Exception:
                traceback.print_exc()
                all_warnings.append("Error running statistical tests")
        else:
            all_warnings.append("Need at least 2 groups for statistical testing")

    # Distance matrix as lists
    dm_data = {
        "labels": list(dm.ids),
        "matrix": dm.data.tolist(),
    }

    # Ordination warnings
    if ord_method == "nmds" and ordination.get("stress") and ordination["stress"] > 0.2:
        all_warnings.append(
            f"NMDS stress = {ordination['stress']} (> 0.2) — results may not reliably represent distances"
        )
    if ord_method == "tsne":
        all_warnings.append(
            "t-SNE is for exploration only — do not use for hypothesis testing"
        )

    return jsonify({
        "samples": samples,
        "distance_metric": metric,
        "ordination_method": ord_method,
        "ordination": ordination,
        "distance_matrix": dm_data,
        "metadata": metadata_dict if metadata_dict else None,
        "metadata_columns": metadata_columns,
        "statistics": statistics,
        "rarefaction_depth": rarefaction_depth,
        "rarefaction_seed": seed,
        "sample_depths": sample_depths,
        "num_samples": len(samples),
        "num_features": int(abundance_df.shape[0]),
        "warnings": all_warnings,
    })


@diversity_bp.route("/api/diversity/mantel", methods=["POST"])
def api_mantel():
    """Run Mantel test between two distance matrices."""
    import io
    import pandas as pd
    from skbio import DistanceMatrix

    if "distance_matrix_1" not in request.files or "distance_matrix_2" not in request.files:
        return jsonify({"error": "Two distance matrix CSV files required"}), 400

    permutations = int(request.form.get("permutations", 999))

    try:
        def parse_dm(file_obj):
            text = file_obj.read().decode("utf-8")
            df = pd.read_csv(io.StringIO(text), index_col=0)
            return DistanceMatrix(df.values, ids=list(df.index.astype(str)))

        dm1 = parse_dm(request.files["distance_matrix_1"])
        dm2 = parse_dm(request.files["distance_matrix_2"])

        if set(dm1.ids) != set(dm2.ids):
            return jsonify({"error": "Distance matrices must have identical sample labels"}), 400

        result = run_mantel(dm1, dm2, permutations=permutations)
        return jsonify(result)

    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    except Exception:
        traceback.print_exc()
        return jsonify({"error": "Error running Mantel test"}), 500
