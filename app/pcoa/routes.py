"""PCoA routes."""

import traceback

from flask import jsonify, render_template, request

from . import pcoa_bp
from .utils import DISTANCE_METRICS, NORMALISATION_METHODS, run_pipeline

import os

BASE_URL = os.environ.get("BASE_URL", "http://localhost:5590").rstrip("/")


@pcoa_bp.route("/pcoa")
def pcoa_page():
    return render_template(
        "pcoa.html",
        distance_metrics=DISTANCE_METRICS,
        normalisation_methods=NORMALISATION_METHODS,
        base_url=BASE_URL,
    )


@pcoa_bp.route("/api/pcoa", methods=["POST"])
def api_pcoa():
    """Run PCoA analysis on an abundance table."""
    # Abundance table (required)
    if "abundance_file" not in request.files or not request.files["abundance_file"].filename:
        return jsonify({"error": "No abundance table uploaded"}), 400

    abundance_bytes = request.files["abundance_file"].read()

    # Metadata (optional)
    metadata_bytes = None
    if "metadata_file" in request.files and request.files["metadata_file"].filename:
        metadata_bytes = request.files["metadata_file"].read()

    metric = request.form.get("distance_metric", "braycurtis")
    norm_method = request.form.get("normalisation", "none")

    try:
        result = run_pipeline(abundance_bytes, metadata_bytes, metric, norm_method)
        return jsonify(result)
    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    except Exception:
        traceback.print_exc()
        return jsonify({"error": "Unexpected server error"}), 500
