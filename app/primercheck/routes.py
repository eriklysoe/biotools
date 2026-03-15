"""PrimerCheck routes."""

import os
import traceback
import threading

from flask import jsonify, render_template, request

from . import primercheck_bp
from .utils_shared import validate_sequence, convert_oligo_type, parse_fasta_text
from .utils_analyze import analyze_oligo, HAS_PRIMER3
from .utils_structure import calc_hairpin, calc_homodimer, calc_heterodimer, generate_dimer_ascii
from .utils_mismatch import calc_tm_mismatch_heatmap
from .utils_check import run_check
from .utils_design import design_primers
from .utils_blast import run_blast

BASE_URL = os.environ.get("BASE_URL", "http://localhost:5590").rstrip("/")


def _get_salt_params(form):
    """Extract salt parameters from form data."""
    return {
        "na_mM": float(form.get("na_conc", 50)),
        "mg_mM": float(form.get("mg_conc", 0)),
        "dntp_mM": float(form.get("dntp_conc", 0)),
        "oligo_nM": float(form.get("oligo_conc", 250)),
    }


@primercheck_bp.route("/primercheck")
def primercheck_page():
    return render_template(
        "primercheck.html",
        base_url=BASE_URL,
        has_primer3=HAS_PRIMER3,
    )


@primercheck_bp.route("/api/primercheck/analyze", methods=["POST"])
def api_analyze():
    """Full oligo property analysis."""
    data = request.get_json(force=True)
    seq = data.get("sequence", "").strip()
    oligo_type = data.get("oligo_type", "DNA")
    target_type = data.get("target_type", "DNA")

    cleaned, errors = validate_sequence(seq, oligo_type)
    if errors:
        return jsonify({"error": errors[0]}), 400

    salt = _get_salt_params(data)

    try:
        result = analyze_oligo(cleaned, oligo_type=oligo_type,
                               target_type=target_type, **salt)
        return jsonify(result)
    except RuntimeError as e:
        return jsonify({"error": str(e)}), 400
    except Exception:
        traceback.print_exc()
        return jsonify({"error": "Unexpected server error"}), 500


@primercheck_bp.route("/api/primercheck/hairpin", methods=["POST"])
def api_hairpin():
    """Hairpin analysis."""
    data = request.get_json(force=True)
    seq = data.get("sequence", "").strip()
    oligo_type = data.get("oligo_type", "DNA")

    cleaned, errors = validate_sequence(seq, oligo_type)
    if errors:
        return jsonify({"error": errors[0]}), 400

    if len(cleaned) > 200:
        warning = "Sequence > 200 nt — structure prediction may be unreliable"
    else:
        warning = None

    salt = _get_salt_params(data)

    try:
        result = calc_hairpin(cleaned, **salt)
        if warning:
            result["length_warning"] = warning
        return jsonify(result)
    except RuntimeError as e:
        return jsonify({"error": str(e)}), 400
    except Exception:
        traceback.print_exc()
        return jsonify({"error": "Unexpected server error"}), 500


@primercheck_bp.route("/api/primercheck/selfdimer", methods=["POST"])
def api_selfdimer():
    """Self-dimer (homodimer) analysis."""
    data = request.get_json(force=True)
    seq = data.get("sequence", "").strip()
    oligo_type = data.get("oligo_type", "DNA")

    cleaned, errors = validate_sequence(seq, oligo_type)
    if errors:
        return jsonify({"error": errors[0]}), 400

    salt = _get_salt_params(data)

    try:
        result = calc_homodimer(cleaned, **salt)
        result["alignment"] = generate_dimer_ascii(
            cleaned.replace("U", "T"), cleaned.replace("U", "T"),
            is_homodimer=True,
        )
        return jsonify(result)
    except RuntimeError as e:
        return jsonify({"error": str(e)}), 400
    except Exception:
        traceback.print_exc()
        return jsonify({"error": "Unexpected server error"}), 500


@primercheck_bp.route("/api/primercheck/heterodimer", methods=["POST"])
def api_heterodimer():
    """Hetero-dimer analysis between two oligos."""
    data = request.get_json(force=True)
    seq1 = data.get("sequence", "").strip()
    seq2 = data.get("sequence2", "").strip()
    oligo_type = data.get("oligo_type", "DNA")

    cleaned1, errors1 = validate_sequence(seq1, oligo_type)
    if errors1:
        return jsonify({"error": f"Sequence 1: {errors1[0]}"}), 400

    cleaned2, errors2 = validate_sequence(seq2, oligo_type)
    if errors2:
        return jsonify({"error": f"Sequence 2: {errors2[0]}"}), 400

    salt = _get_salt_params(data)

    try:
        result = calc_heterodimer(cleaned1, cleaned2, **salt)
        result["alignment"] = generate_dimer_ascii(
            cleaned1.replace("U", "T"), cleaned2.replace("U", "T"),
        )
        return jsonify(result)
    except RuntimeError as e:
        return jsonify({"error": str(e)}), 400
    except Exception:
        traceback.print_exc()
        return jsonify({"error": "Unexpected server error"}), 500


@primercheck_bp.route("/api/primercheck/tmmismatch", methods=["POST"])
def api_tmmismatch():
    """Tm mismatch heatmap."""
    data = request.get_json(force=True)
    seq = data.get("sequence", "").strip()
    oligo_type = data.get("oligo_type", "DNA")

    cleaned, errors = validate_sequence(seq, oligo_type)
    if errors:
        return jsonify({"error": errors[0]}), 400

    salt = _get_salt_params(data)

    try:
        result = calc_tm_mismatch_heatmap(cleaned, **salt)
        return jsonify(result)
    except RuntimeError as e:
        return jsonify({"error": str(e)}), 400
    except Exception:
        traceback.print_exc()
        return jsonify({"error": "Unexpected server error"}), 500


@primercheck_bp.route("/api/primercheck/check", methods=["POST"])
def api_check():
    """In-silico PCR: check primers against template."""
    data = request.get_json(force=True)
    fwd = data.get("forward_primer", "").strip()
    rev = data.get("reverse_primer", "").strip()
    template = data.get("template", "").strip()
    max_mm = int(data.get("max_mismatches", 2))
    oligo_type = data.get("oligo_type", "DNA")

    # Validate primers
    fwd_clean, fwd_err = validate_sequence(fwd, oligo_type)
    if fwd_err:
        return jsonify({"error": f"Forward primer: {fwd_err[0]}"}), 400
    rev_clean, rev_err = validate_sequence(rev, oligo_type)
    if rev_err:
        return jsonify({"error": f"Reverse primer: {rev_err[0]}"}), 400

    # Parse template (could be FASTA or plain)
    if not template:
        return jsonify({"error": "No template sequence provided"}), 400

    # If FASTA, take first sequence
    sequences = parse_fasta_text(template)
    if not sequences:
        return jsonify({"error": "Could not parse template sequence"}), 400

    if len(sequences) > 10:
        template_warning = "Template has >10 sequences; using first 10"
        sequences = sequences[:10]
    else:
        template_warning = None

    salt = _get_salt_params(data)

    try:
        # Run check against first template sequence
        template_seq = sequences[0][1]
        result = run_check(
            fwd_clean, rev_clean, template_seq,
            max_mismatches=max_mm, oligo_type=oligo_type, **salt,
        )
        result["template_name"] = sequences[0][0]
        if template_warning:
            result["warnings"].append(template_warning)
        return jsonify(result)
    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    except RuntimeError as e:
        return jsonify({"error": str(e)}), 400
    except Exception:
        traceback.print_exc()
        return jsonify({"error": "Unexpected server error"}), 500


@primercheck_bp.route("/api/primercheck/design", methods=["POST"])
def api_design():
    """Design primer pairs using Primer3."""
    data = request.get_json(force=True)
    template = data.get("template", "").strip()

    if not template:
        return jsonify({"error": "No template sequence provided"}), 400

    # Parse template
    sequences = parse_fasta_text(template)
    if not sequences:
        return jsonify({"error": "Could not parse template sequence"}), 400

    template_seq = sequences[0][1]
    salt = _get_salt_params(data)

    # Parse excluded regions
    excluded = None
    excl_text = data.get("excluded_regions", "").strip()
    if excl_text:
        try:
            pairs = excl_text.split(",")
            excluded = []
            for i in range(0, len(pairs), 2):
                excluded.append([int(pairs[i].strip()), int(pairs[i+1].strip())])
        except (ValueError, IndexError):
            return jsonify({"error": "Invalid excluded regions format. Use: start,length,start,length,..."}), 400

    # Target region
    target_start = data.get("target_start")
    target_length = data.get("target_length")
    if target_start is not None and target_length is not None:
        try:
            target_start = int(target_start)
            target_length = int(target_length)
        except (ValueError, TypeError):
            target_start = None
            target_length = None

    try:
        result = design_primers(
            template_seq,
            target_start=target_start,
            target_length=target_length,
            excluded_regions=excluded,
            min_size=int(data.get("min_size", 100)),
            max_size=int(data.get("max_size", 500)),
            min_length=int(data.get("min_length", 18)),
            opt_length=int(data.get("opt_length", 20)),
            max_length=int(data.get("max_length", 25)),
            min_tm=float(data.get("min_tm", 57)),
            opt_tm=float(data.get("opt_tm", 60)),
            max_tm=float(data.get("max_tm", 63)),
            min_gc=float(data.get("min_gc", 40)),
            max_gc=float(data.get("max_gc", 60)),
            max_tm_diff=float(data.get("max_tm_diff", 3)),
            num_pairs=int(data.get("num_pairs", 5)),
            **salt,
        )
        result["template_name"] = sequences[0][0]
        return jsonify(result)
    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    except RuntimeError as e:
        return jsonify({"error": str(e)}), 400
    except Exception:
        traceback.print_exc()
        return jsonify({"error": "Unexpected server error"}), 500


@primercheck_bp.route("/api/primercheck/blast", methods=["POST"])
def api_blast():
    """BLAST oligo against NCBI."""
    data = request.get_json(force=True)
    seq = data.get("sequence", "").strip()
    oligo_type = data.get("oligo_type", "DNA")

    cleaned, errors = validate_sequence(seq, oligo_type)
    if errors:
        return jsonify({"error": errors[0]}), 400

    try:
        result = run_blast(cleaned)
        return jsonify(result)
    except TimeoutError as e:
        return jsonify({"error": str(e), "timeout": True}), 408
    except requests.exceptions.HTTPError as e:
        if e.response is not None and e.response.status_code == 429:
            return jsonify({
                "error": "NCBI rate limit reached. Please retry after 60 seconds.",
                "rate_limited": True,
            }), 429
        return jsonify({"error": f"BLAST HTTP error: {e}"}), 502
    except Exception:
        traceback.print_exc()
        return jsonify({"error": "Unexpected server error during BLAST"}), 500
