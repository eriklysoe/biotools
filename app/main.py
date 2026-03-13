"""
BioTools – Combined bioinformatics web tools.
Includes SeqConvert (sequence format converter) and VennApp (Venn diagram tool).
"""

import os
import uuid
import tempfile
import traceback
from pathlib import Path
from functools import wraps

from flask import (
    Flask, request, render_template, send_file,
    jsonify, abort, Response,
)
from werkzeug.middleware.proxy_fix import ProxyFix
from Bio import SeqIO

from .formats import INPUT_FORMATS, OUTPUT_FORMATS, MOLECULE_TYPES, FORMAT_INFO

app = Flask(
    __name__,
    template_folder=os.path.join(os.path.dirname(__file__), "..", "templates"),
    static_folder=os.path.join(os.path.dirname(__file__), "..", "static"),
)

# Trust X-Forwarded-* headers from reverse proxy
app.wsgi_app = ProxyFix(app.wsgi_app, x_for=1, x_proto=1, x_host=1, x_prefix=1)

MAX_UPLOAD_MB = int(os.environ.get("MAX_UPLOAD_MB", 50))
BASE_URL = os.environ.get("BASE_URL", "http://localhost:5590").rstrip("/")
app.config["MAX_CONTENT_LENGTH"] = MAX_UPLOAD_MB * 1024 * 1024

TMP_DIR = Path(os.environ.get("TMP_DIR", "/app/tmp"))
TMP_DIR.mkdir(parents=True, exist_ok=True)

ADMIN_USER = os.environ.get("ADMIN_USER", "")
ADMIN_PASS = os.environ.get("ADMIN_PASS", "")


def require_auth(f):
    """Decorator: enforce HTTP Basic Auth if ADMIN_USER and ADMIN_PASS are set."""
    @wraps(f)
    def decorated(*args, **kwargs):
        if not ADMIN_USER or not ADMIN_PASS:
            return f(*args, **kwargs)
        auth = request.authorization
        if not auth or auth.username != ADMIN_USER or auth.password != ADMIN_PASS:
            return Response(
                "Login required.", 401,
                {"WWW-Authenticate": 'Basic realm="BioTools"'},
            )
        return f(*args, **kwargs)
    return decorated


# ── SeqConvert logic ─────────────────────────────────────────────────

FORMAT_EXTENSIONS = {
    "fasta": ".fasta", "fasta-2line": ".fasta",
    "fastq": ".fastq", "fastq-sanger": ".fastq",
    "fastq-solexa": ".fastq", "fastq-illumina": ".fastq",
    "genbank": ".gb", "gb": ".gb",
    "embl": ".embl", "imgt": ".embl",
    "clustal": ".aln", "nexus": ".nex",
    "phylip": ".phy", "stockholm": ".sto",
    "seqxml": ".xml", "phd": ".phd",
    "pir": ".pir", "tab": ".tsv",
    "qual": ".qual", "sff": ".sff", "xdna": ".xdna",
}


def _annotate_molecule_type(records, mol_type):
    """Yield records with molecule_type annotation set."""
    for record in records:
        if mol_type:
            record.annotations["molecule_type"] = mol_type
        yield record


def convert_sequences(input_path, in_fmt, out_fmt, mol_type_key="none"):
    """Convert a sequence file and return (output_path, record_count)."""
    if in_fmt not in INPUT_FORMATS:
        raise ValueError(f"Unsupported input format: {in_fmt}")
    if out_fmt not in OUTPUT_FORMATS:
        raise ValueError(f"Unsupported output format: {out_fmt}")

    mol_type = MOLECULE_TYPES.get(mol_type_key)
    ext = FORMAT_EXTENSIONS.get(out_fmt, ".txt")
    out_name = f"{uuid.uuid4().hex}{ext}"
    output_path = TMP_DIR / out_name

    try:
        if mol_type:
            records = SeqIO.parse(str(input_path), in_fmt)
            annotated = _annotate_molecule_type(records, mol_type)
            count = SeqIO.write(annotated, str(output_path), out_fmt)
        else:
            count = SeqIO.convert(str(input_path), in_fmt, str(output_path), out_fmt)
    except Exception as e:
        output_path.unlink(missing_ok=True)
        raise RuntimeError(f"Conversion failed: {e}") from e

    if count == 0:
        output_path.unlink(missing_ok=True)
        raise RuntimeError(
            "No records were converted. Check that the input file matches "
            f"the selected format ({in_fmt})."
        )

    return output_path, count


# ── Routes: SeqConvert ───────────────────────────────────────────────

@app.route("/")
@app.route("/seqconvert")
@require_auth
def index():
    return render_template(
        "seqconvert.html",
        input_formats=INPUT_FORMATS,
        output_formats=OUTPUT_FORMATS,
        molecule_types=list(MOLECULE_TYPES.keys()),
        format_info=FORMAT_INFO,
        base_url=BASE_URL,
    )


@app.route("/api/convert", methods=["POST"])
@require_auth
def api_convert():
    """Convert an uploaded sequence file."""
    if "file" not in request.files:
        return jsonify({"error": "No file uploaded"}), 400

    f = request.files["file"]
    if f.filename == "":
        return jsonify({"error": "Empty filename"}), 400

    in_fmt = request.form.get("input_format", "fasta")
    out_fmt = request.form.get("output_format", "fasta")
    mol_type = request.form.get("molecule_type", "none")

    suffix = Path(f.filename).suffix or ".seq"
    tmp_in = tempfile.NamedTemporaryFile(
        dir=TMP_DIR, suffix=suffix, delete=False
    )
    try:
        f.save(tmp_in)
        tmp_in.close()

        output_path, count = convert_sequences(tmp_in.name, in_fmt, out_fmt, mol_type)

        stem = Path(f.filename).stem
        ext = FORMAT_EXTENSIONS.get(out_fmt, ".txt")
        download_name = f"{stem}_converted{ext}"

        response = send_file(
            str(output_path),
            as_attachment=True,
            download_name=download_name,
        )

        @response.call_on_close
        def _cleanup():
            Path(output_path).unlink(missing_ok=True)

        return response

    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    except RuntimeError as e:
        return jsonify({"error": str(e)}), 422
    except Exception:
        traceback.print_exc()
        return jsonify({"error": "Unexpected server error"}), 500
    finally:
        Path(tmp_in.name).unlink(missing_ok=True)


@app.route("/api/formats", methods=["GET"])
@require_auth
def api_formats():
    """Return available formats as JSON."""
    return jsonify({
        "input_formats": INPUT_FORMATS,
        "output_formats": OUTPUT_FORMATS,
        "molecule_types": list(MOLECULE_TYPES.keys()),
        "format_info": FORMAT_INFO,
    })


# ── Routes: VennApp ──────────────────────────────────────────────────

@app.route("/venn")
@require_auth
def venn():
    return render_template("venn.html")


@app.route("/api/compare", methods=["POST"])
@require_auth
def api_compare():
    """Compare ID sets and return intersection regions."""
    data = request.get_json(force=True)
    sets_raw = data.get("sets", {})

    if len(sets_raw) < 2 or len(sets_raw) > 4:
        return jsonify({"error": "Provide between 2 and 4 sets"}), 400

    named_sets = {name: set(ids) for name, ids in sets_raw.items()}
    names = list(named_sets.keys())
    n = len(names)

    results = []
    for mask in range(1, 1 << n):
        included = [names[i] for i in range(n) if mask & (1 << i)]
        excluded = [names[i] for i in range(n) if not (mask & (1 << i))]

        region = set.intersection(*(named_sets[nm] for nm in included))
        for nm in excluded:
            region = region - named_sets[nm]

        results.append({
            "included": included,
            "excluded": excluded,
            "ids": sorted(region),
            "count": len(region),
        })

    return jsonify({"regions": results, "set_names": names})


# ── Health check ─────────────────────────────────────────────────────

@app.route("/health")
def health():
    return jsonify({"status": "ok"})
