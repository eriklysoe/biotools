"""PhyloTree routes."""

import shutil
import subprocess
import tempfile
import traceback
import uuid
from pathlib import Path

from flask import jsonify, render_template, request, Response

from . import phylotree_bp
from .utils import EXPORT_FORMATS, get_available_methods, newick_to_format, run_pipeline

import os

TMP_DIR = Path(os.environ.get("TMP_DIR", "/app/tmp"))
BASE_URL = os.environ.get("BASE_URL", "http://localhost:5590").rstrip("/")

MAX_SEQUENCES_WARN = 200


@phylotree_bp.route("/phylotree")
def phylotree_page():
    available = get_available_methods()
    return render_template(
        "phylotree.html",
        available=available,
        base_url=BASE_URL,
    )


@phylotree_bp.route("/api/phylotree", methods=["POST"])
def api_phylotree():
    """Build a phylogenetic tree from FASTA input."""
    # Accept file upload OR pasted text
    fasta_text = None
    original_name = "sequences.fasta"

    if "file" in request.files and request.files["file"].filename:
        f = request.files["file"]
        original_name = f.filename
        fasta_text = f.read().decode("utf-8", errors="replace")
    elif request.form.get("sequences", "").strip():
        fasta_text = request.form["sequences"].strip()
    else:
        return jsonify({"error": "No FASTA file or sequences provided"}), 400

    align_method = request.form.get("alignment_method", "muscle")
    tree_method = request.form.get("tree_method", "nj")

    # Validate input looks like FASTA
    if not fasta_text.startswith(">"):
        return jsonify({
            "error": "Invalid FASTA format: sequences must start with '>'"
        }), 400

    # Create a unique working directory
    job_id = uuid.uuid4().hex
    work_dir = TMP_DIR / f"phylo_{job_id}"
    work_dir.mkdir(parents=True, exist_ok=True)

    input_path = work_dir / "input.fasta"

    try:
        input_path.write_text(fasta_text)

        # Quick sequence count check
        seq_count = fasta_text.count("\n>") + 1
        if seq_count > MAX_SEQUENCES_WARN and tree_method == "ml":
            return jsonify({
                "error": (
                    f"Dataset has ~{seq_count} sequences. "
                    "ML analysis on large datasets may be very slow. "
                    "Consider using NJ/UPGMA or running IQ-TREE locally."
                ),
                "warning": True,
            }), 422

        result = run_pipeline(input_path, align_method, tree_method, work_dir)
        return jsonify(result)

    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    except subprocess.TimeoutExpired:
        return jsonify({
            "error": "Analysis timed out (5 min limit). Try fewer sequences or a faster method."
        }), 422
    except RuntimeError as e:
        return jsonify({"error": str(e)}), 422
    except Exception:
        traceback.print_exc()
        return jsonify({"error": "Unexpected server error"}), 500
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


@phylotree_bp.route("/api/phylotree/export", methods=["POST"])
def api_export():
    """Convert a Newick string to another tree file format."""
    data = request.get_json(force=True)
    newick = data.get("newick", "").strip()
    fmt = data.get("format", "newick")

    if not newick:
        return jsonify({"error": "No Newick string provided"}), 400

    if fmt == "newick":
        content = newick
    else:
        try:
            content = newick_to_format(newick, fmt)
        except ValueError as e:
            return jsonify({"error": str(e)}), 400
        except Exception:
            traceback.print_exc()
            return jsonify({"error": "Failed to convert tree format"}), 500

    info = EXPORT_FORMATS.get(fmt, EXPORT_FORMATS["newick"])
    return Response(
        content,
        mimetype=info["mime"],
        headers={
            "Content-Disposition": f'attachment; filename="tree{info["ext"]}"'
        },
    )


@phylotree_bp.route("/api/phylotree/methods")
def api_methods():
    """Return available alignment and tree methods."""
    return jsonify(get_available_methods())
